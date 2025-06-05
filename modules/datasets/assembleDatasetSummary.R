assembleDatasetSummary <- function(input){
  
  datasetName <- metadata[input$metadataTable_rows_selected,]$Folder
  
  datasetSummary <- datasetDir %>%
    list.files(full.names = TRUE, recursive = TRUE) %>%
    {.[map_lgl(., function(x){grepl("meta_data.md", x)})]} %>%
    {.[str_detect(., datasetName)]} %>%
    pandoc_convert(to = "html") %>%
    paste() %>%
    HTML

  return(datasetSummary)
  
}