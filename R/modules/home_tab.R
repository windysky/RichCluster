
homeTabUI <- function(id, tabName) {
  ns <- shiny::NS(id)
  shiny::tabItem(tabName = tabName,
          shiny::h1("RichCluster Shiny App"),
          shiny::p("Developed by Sarah Hong")
  )
}


homeTabServer <- function(id) {

  shiny::moduleServer(id, function(input, output, session) {


  })

}
