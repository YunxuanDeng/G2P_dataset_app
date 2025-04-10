source("global.R")
source("RSFunctions.R")

ui <- navbarPage(id = "tabs", collapsible = TRUE, title = "Genome-to-Phenome Datasets App",
  modulePanel("Datasets", value = "datasets"),
  modulePanel("Add dataset", value = "addDataset"),
  modulePanel("About", value = "about"),
  reloadWarning,
  tags$head(tags$link(rel="shortcut icon", href="https://storage.googleapis.com/kaggle-organizations/4421/thumbnail.png?t=2025-01-17-04-03-19"))
)
  
server <- function(input, output, session) {
  
  sapply(list.files(path = "modules"), function(module){
    get(module)$server(input, output, session)
  })
  
}

shinyApp(ui = ui, server = server)