assembleDatasetCitation <- function(datasets){
  
  datasets <- metadata[datasets,]$Folder
  datasetDir <- paste0("GPDatasets-main/Datasets/", datasets, "/output")
    # Need to modify to be automatically set to "testing-main/Datasets/{Folder}/{ouput}"
  
  Citation <- datasetDir %>%
    list.files(full.names = TRUE, recursive = TRUE) %>%
    {.[map_lgl(., function(x){grepl("citation.bib", x)})]} %>%
    {.[str_detect(., paste(datasets, collapse = "|"))]} %>%
    map(., readLines) %>%
    map(., ~paste(.x, collapse = "\n")) %>%
    paste(., collapse = "\n")
  
  return(Citation)
  
}