
visualizeTabUI <- function(id, tabName) {
  ns <- shiny::NS(id)
  shiny::tabItem(tabName = tabName,
          shiny::h1("Welcome to RichStudio v0.1.5!"),
          shiny::p("Upload from differential expression geneset (DEG) or enrichment result. Supports kappa clustering and multiple visualizations.")
  )
}


visualizeTabServer <- function(id) {

  shiny::moduleServer(id, function(input, output, session) {


  })

}
