assembleCurationCode <- function(datasets){
  
  datasets <- metadata[datasets,]$Folder
  
  CurationCode <- datasetDir %>%
    list.files(full.names = TRUE, recursive = TRUE) %>%
    {.[map_lgl(., function(x){grepl("curate_data_code.R", x)})]} %>%
    {.[str_detect(., paste(datasets, collapse = "|"))]} %>%
    map(., readLines) %>%
    map(., ~paste(.x, collapse = "\n")) %>%
    paste(., collapse = "\n")
  
  return(CurationCode)
  
}