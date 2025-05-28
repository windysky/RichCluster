
clusterTabUI <- function(id, tabName) {
  ns <- shiny::NS(id)
  shiny::tabItem(tabName = tabName,
    shiny::h1("Cluster Enrichment Results"),
    shiny::fluidRow(
      shiny::column(width=4,
        shinydashboard::box(
          title = "Upload Enrichment Results",
          status = "primary",
          solidHeader = TRUE,
          collapsible = TRUE,
          width = NULL,
          shiny::fileInput(ns("rich_upload"), 'Select files', multiple=TRUE, accept=c('.csv', '.tsv', '.xls', '.txt')),
          shiny::actionButton(ns("submit"), "Submit"),
          shiny::hr(),
          shiny::h4('Load demo'),
          shiny::p("Loads sample files: 'go1', 'go2'"),
          shiny::actionButton(ns("rich_load"), "Load")
        )
      ),
      shiny::column(width=8,
        shinydashboard::tabBox(width=NULL,
          shiny::tabPanel("Distance Matrix",
            shiny::h4("Select enrichment results"),
            shiny::selectInput(ns('select_rs_for_distmatrix'), "Select enrichment results to compute distances", choices=NULL, multiple=TRUE),
            shiny::selectInput(ns('dist_metric'), "Select distance metric", choices=c("kappa", "jaccard")),
            shiny::selectInput(ns('membership_strategy'), "Select membership strategy", choices=c("DAVID", "Ward")),
            shiny::actionButton(ns('compute_distances'), 'Compute Distances'),
            shiny::br(),
            shiny::br(),
            DT::dataTableOutput(ns("distance_matrix"))
          ),
          shiny::tabPanel("Filter Clusters",
            shiny::numericInput(ns('min_terms'), "Minimum number of terms", value=5),
            shiny::actionButton(ns('filter_clusters'), 'Filter Clusters'),
            DT::dataTableOutput(ns("merged_seeds"))
          ),
          shiny::tabPanel("Cluster Correlation",
            shiny::selectInput(ns('select_cluster'), "Select cluster", choices=NULL),
            plotly::plotlyOutput(ns("cluster_correlation"))
          ),
          shiny::tabPanel("Full Heatmap",
            plotly::plotlyOutput(ns("full_heatmap"))
          ),
        )
      )
    )
  )
}


clusterTabServer <- function(id, richnames, richsets) {

  shiny::moduleServer(id, function(input, output, session) {

    select_richnames <- shiny::reactive(richnames$labels)
    merged_richsets <- shiny::reactiveVal(NULL)

    distance_matrix <- shiny::reactiveVal(NULL)
    merged_seeds <- shiny::reactiveVal(NULL)

    final_clusters <- shiny::reactiveVal(NULL)
    full_clusterdf <- shiny::reactiveVal(NULL)

    shiny::observe({
      shiny::updateSelectInput(session=shiny::getDefaultReactiveDomain(), 'select_rs_for_distmatrix', choices=select_richnames())
    })

    # <!----- UPLOAD LOGIC -----!>
    shiny::observeEvent(input$rich_upload, {
      shiny::req(input$rich_upload) # Make sure file uploaded

      for (i in seq_along(input$rich_upload$name)) {
        richname <- input$rich_upload$name[i]

        file_extension <- tools::file_ext(input$rich_upload$name[i]) # tools::file_ext already correct
        file_path <- input$rich_upload$datapath[i]

        # try to read file as csv
        csv_ncol <- tryCatch({
          csvdf <- read.csv(path)
          ncol(csvdf)
        }, error = function(err) {
          0
        })
        # try to read file as tsv
        tsv_ncol <- tryCatch({
          tsvdf <- read.delim(path)
          ncol(tsvdf)
        }, error = function(err) {
          0
        })
        # decide which df to store
        if (tsv_ncol == 0 || csv_ncol > tsv_ncol) {
          rich_df <- read.csv(path)
        } else {
          rich_df <- read.delim(path)
        }

        #u_rrdfs[[lab]] <- select_required_columns(df) # set u_rrdfs
        richsets[[richname]] <- rich_df # set u_rrdfs
        richnames$labels <- c(richnames$labels, richname) # set richnames
      }

    })

    observeEvent(input$rich_load, {
      # go1 and go2
      for (i in 1:2) {
        file_name <- paste0('go', i)
        file_path <- system.file('extdata', paste0(file_name, '.txt'), package='RichCluster')
        rich_df <- read.delim(file_path)

        richsets[[file_name]] <- rich_df
        richnames$labels <- c(richnames$labels, file_name)
      }
      shiny::showNotification("Loaded demo files")
    })

    shiny::observeEvent(input$compute_distances, {
      shiny::req(input$select_rs_for_distmatrix)
      shiny::req(input$dist_metric)

      shiny::withProgress(message="Computing distances...", value=0, {
        selected_richnames <- input$select_rs_for_distmatrix
        rich_dfs <- lapply(selected_richnames, function(richname) {
          richsets[[richname]]
        })

        merged_richsets <- RichCluster::merge_richsets(rich_dfs)
        term_vec <- merged_richsets$Term
        geneID_vec <- merged_richsets$GeneID
        padj_vec <- merged_richsets$Padj

        shiny::incProgress(0.3, message=NULL)

        cluster_results <- tryCatch(
          RichCluster::RichCluster(
            "kappa",
            0.5,
            "DAVID",
            0.5,
            term_vec,
            geneID_vec,
            padj_vec
          ),
          error = function(err) {
            shiny::showNotification(err$message, type='error')
            shiny::incProgress(0.7)
            return(NULL)
          }
        )
        shiny::incProgress(0.6, message=NULL)

        if (!is.null(cluster_results)) {
          distance_matrix(cluster_results$DistanceMatrix)
          merged_seeds(cluster_results$MergedSeeds)

          distance_matrix(as.data.frame(distance_matrix))
          output$distance_matrix <- DT::renderDataTable(distance_matrix())

          output$merged_seeds <- DT::renderDataTable(merged_seeds())

          final_clusters(RichCluster::filter_merged_seeds(merged_seeds(), 5))
          full_clusterdf(RichCluster::make_full_clusterdf(final_clusters(), merged_richsets()))

        }
        shiny::incProgress(0.1, message=NULL)
        shiny::showNotification("Done!")
      })
    })

    shiny::observeEvent(input$filter_clusters, {
      shiny::req(input$min_terms)
      shiny::req(merged_seeds())

      final_clusters(RichCluster::filter_merged_seeds(merged_seeds(), input$min_terms))
      full_clusterdf(RichCluster::make_full_clusterdf(final_clusters(), merged_richsets()))

      # output$merged_seeds <- DT::renderDataTable(final_clusters)
      # output$full_heatmap <- plotly::renderPlotly({
      #   RichCluster::all_clusters_hmap(final_clusters(), "Padj")
      # })
    })

    # observe({
    #   if (!is.null(final_clusters())) {
    #     output$full_heatmap <- plotly::renderPlotly({
    #       RichCluster::all_clusters_hmap(final_clusters(), "Padj")
    #     })
    #   }
    # })

    table_merged_seeds <- shiny::reactive({
      shiny::req(!is.null(full_clusterdf()))
      return(full_clusterdf())
    })
    output$merged_seeds <- DT::renderDT({
      table_merged_seeds()
    })

    plot_all_clusters_hmap <- shiny::reactive({
      shiny::req(!is.null(full_clusterdf()))
      hmap <- RichCluster::all_clusters_hmap(full_clusterdf(), "Padj")
      return(hmap)
    })
    output$full_heatmap = plotly::renderPlotly({ # plotly::renderPlotly already correct
      hmap <- RichCluster::all_clusters_hmap(full_clusterdf(), "Padj")
    })

  })

}
