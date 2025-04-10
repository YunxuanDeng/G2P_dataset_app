assembleAnalysesCode <- function(datasets){
  
  AnalysesCodeDir <- "GPDatasets-main/Analyses"
  datasets <- metadata[datasets,]$Folder
  
  AnalysesCode <- AnalysesCodeDir %>%
    list.files(full.names = TRUE, recursive = TRUE) %>%
    map_chr(function(file_path) {
      # Extract method name from the file name
      method_name <- str_extract(file_path, "(?<=Analyses/).*(?=\\.R)")
      
      # Read the file content
      file_content <- suppressWarnings(readLines(file_path))
      
      # Add the method name as a comment at the beginning
      commented_content <- c(paste0("# Method: ", method_name), file_content)
      
      # Collapse the content into a single string for return
      paste(commented_content, collapse = "\n")
    }) %>%
    paste(collapse = "\n\n")
  
  return(AnalysesCode)
}