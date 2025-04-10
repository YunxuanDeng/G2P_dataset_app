assembleMetadata <- function(metadata){

  createLink <- function(val) {
    sprintf(paste0('<a href="', URLdecode(val),'" target="_blank">', val,'</a>'))
  }

  completedDatasets <- datasetDir %>%
    list.files(full.names = TRUE, recursive = TRUE) %>%
    {.[map_lgl(., function(x){grepl("curate_data_code.R", x)})]} %>%
    {.[file.info(.)$size > 1]}
  
  metadata %<>%
    rowwise() %>%
    mutate(
      `In Kaggle?` = any(grepl(Folder, completedDatasets)),
      Kaggle = ifelse(`In Kaggle?`, createLink(paste0("https://www.kaggle.com/datasets/ag2p-disc/", gsub("_", "-", tolower(Folder)))), "")
    ) %>%
    select(Title, `N Geno` = N_Geno, `N Markers` = N_markers, `Has Phenotype?` = Phenotype, Kaggle) %>%
    rowwise() %>%
    mutate(
      `N Geno` = as.integer(`N Geno`),
      `N Markers` = as.integer(`N Markers`),
      #DOI = createLink(DOI)
    ) %>%
    ungroup()
  
  return(metadata)
  
}
