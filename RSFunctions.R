require(httr)
require(RCurl)
require(rjson)
library(data.table)
library(DT)
require(rclipboard)
require(tidyverse)
require(magrittr)
library(rjson)
require(jsonlite)
require(curl)
library(readxl)
library(gert)
library(stringdist)
library(base64enc)
require(assignPOP)


# Check if the URL is valid
check_url <- function(url,Meta_Data) {
  response <- tryCatch({GET(url)}, error = function(e) {NULL})

  if (!is.null(response) && (status_code(response) == 200 || status_code(response) == 403) && grepl("https://doi.org/", url)) {
    data_link <- GET(url)$url
    message0 <- "<span style='color: green;'>DOI verified!</span><br>"
    if (!(url %in% Meta_Data$DOI)) {
      message1 <- "<span style='color: black;'>Please verify the data sharing link, modify it if needed and validate it.</span><br>"
      return(list(paste(message0,message1), TRUE, data_link, 0))
    } else{
      result <- Meta_Data[Meta_Data$DOI %in% url, ][,c(1,4,7,13,14)]
      message00 <- "<span style='color: red;'>Data already exists!</span><br>"
      message11 <- paste0("<span style='color: black;'>If you are uploading data that is different from the existing data, please proceed by verifying the data sharing link, modify it if needed and validate it.</span>")
      return(list(paste(message0,message00,message11), TRUE, data_link, 1, result))
    } 
  } else {
    return(list("<span style='color: red;'>Please enter a valid DOI!</span>", FALSE))
  }
} 

check_link <- function(url){
  response <- tryCatch({GET(url)}, error = function(e) {NULL})

  if (!is.null(response) && (status_code(response) == 200 || status_code(response) == 403)) {
    return(list("<span style='color: green;'>Link verified!</span>", TRUE))
  } else{
    return(list("<span style='color: red;'>Please enter a valid link!</span>", FALSE))
  }
}


# Check if the article exists
check_Title <- function(Title,DOI,TF){
  if(!TF){
    if(Title == ""){
      return(list("<span style='color: red;'>Please enter the title of the article!</span>",FALSE))
    } else{
      find_Article(Title,Title)
    }
  } else{
    response <- curl_fetch_memory(DOI, handle = handle_setheaders(new_handle(), accept = "application/x-bibtex"))
    bb <- rawToChar(response$content)
    lines <- unlist(strsplit(bb, "\n"))
    author_line <- grep("author\\s*=\\s*\\{", lines, value = TRUE)
    authors <- sub(".*?\\{(.*?)\\}.*", "\\1", author_line)
    author <- gsub(",","",strsplit(authors, " and ")[[1]][1])
    find_Article(Title,paste(Title,author))
  }
}

# Get article information
find_Article <- function(Title, SearchTerm){
  api_key <- "l8U7jnGFYH6csYDcjy3K31LD3UnBrOYt7HgIW0S7"
  base_url <- "https://api.semanticscholar.org/graph/v1/paper/search/"
  params <- list(
    query = SearchTerm, 
    fields = "title,abstract,externalIds",
    limit = 5
  )
  response <- try(GET(base_url, add_headers(`x-api-key` = api_key), query = params), silent = TRUE)
  result <- fromJSON(content(response, "text", encoding = "UTF-8"))
  if(result$total==0){
    return(list("<span style='color: red;'>The Article couldn't be found!</span>",FALSE)) 
  } else{
    result <- result$data[[1]]
    similarity <- 1 - stringdist(Title, result$title[1], method = "jw")
    if(similarity > 0.8){
      # Get the publication link
      doi_url <- paste0("https://doi.org/",result$externalIds$DOI[1])
      pub_link <- GET(doi_url)$url
      String <- paste("Publication:",pub_link )
      # Get the Abstract
      abstract <- sub("^[\\s:~*]+", "", result$abstract[1])
      message0 <- paste("<span style='color: black;'>Article found at </span><a href='", pub_link, "' target='_blank'>", pub_link, "</a>", sep = "")
      message1 <- "<span style='color: grey; font-size:12px;'>If this is not the desired article, please modify it after previewing the data.</span>"
      return(list(paste(message0,"<br>",message1),TRUE,pub_link,abstract))
    } else{
      return(list("<span style='color: red;'>The Article couldn't be found!</span>",FALSE))
    }
  }
  
}

# ChatGPT
askChatGPT <- function(prompt, chatGPTApiKey, max_tokens = 512){
  response <- httr::POST(
    url = "https://api.openai.com/v1/chat/completions",
    add_headers(Authorization = paste("Bearer", chatGPTApiKey)),
    content_type("application/json"),
    encode = "json",
    body = list(
      model = "gpt-4-turbo",
      messages = list(
        list(role = "user", content = prompt)
      ),
      temperature = 0.7,
      max_tokens = max_tokens
    )
  )
  response <- str_trim(content(response)$choices[[1]]$message$content)
  return(response)
}
askAboutText <- function(Text, task, chatGPTApiKey){
  Text %>%
    paste(task, .) %>%
    substr(1, 10000) %>%
    askChatGPT(., chatGPTApiKey)
}



# Find species
Get_Species <- function(Text){
  chatGPTApiKey <<- "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"
  question <- "What species is the given text is about, Give scientific name and its common name separeted by hypen with their first letters being capital and no other text:"
  askAboutText(Text, question, chatGPTApiKey)
}

# Generate Tags
add_tags <- function(Text){
  chatGPTApiKey <<- "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"
  question <- "Return 5 to 8 tags for the given abstract text and seperate them by semi colons, if there is no text return NA:"
  askAboutText(Text, question, chatGPTApiKey)
}

# Create Meta Data
create_meta_data <- function(DOI,Pub_link,Title,Abstract,Data_Link,Spec,hasPheno,nSamples,nMarkers,Code,Meta_Data){
  # if(Spec==""){
  #   Species <- strsplit(Get_Species(paste(Title,Abstract))," - ")[[1]]
  # } else{
  #   Species <- strsplit(Get_Species(paste(Title,Spec))," - ")[[1]]
  # }
  Species <- c("",paste0(toupper(substring(Spec, 1, 1)), substring(Spec, 2)))
  #Tags <- as.character(add_tags(Abstract))
  Tags <- ""
  folder_name <- generate_folder_name(Title,Species[2],Meta_Data)
  inputs <- c("Unique Name","Data DOI","Article Publication link","Title of the Article","Abtract/Description","Species Scientific Name","Species Common Name",
               "Data Sharing Link","Authorization for Accessing Data","Phenotypic Data","Genetic Map","Pedigree Information","Number of samples",
               "Number of markers","Data Downloading Instructions","Article Tags")
  outputs <- c(folder_name,DOI,Pub_link,Title,Abstract,Species[1],Species[2],Data_Link,FALSE,hasPheno || grepl("Get pheno", Code),
                grepl("Get map", Code),grepl("Get ped", Code),nSamples,nMarkers,NA,Tags)
  return(list(inputs,outputs))
}

# Generate folder name
generate_folder_name <- function(Title,Species_Common_Name,Meta_Data){
  chatGPTApiKey <<- "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"
  question <- "Generate a 3 word concise text for the given text, don't include the name of the species and Capitalize first letter of each word"
  # unique_txt <- askAboutText(Title, question, chatGPTApiKey)
  sentence <- gsub("[^a-zA-Z\\s]", " ", tools::toTitleCase(Title))
  words = strsplit(sentence, " ")[[1]]
  unique_txt <- paste(sort(tail(words[order(nchar(words))],3)), collapse = " ")
  Unique_name <- gsub("[-',:() ]", "", paste(unique_txt, Species_Common_Name))
  Folder_name <- gsub(" ","",paste(sprintf("%05d",nrow(Meta_Data)+1),"_",Unique_name))
  return(Folder_name)
}

# Get Meta data from Github
# Update accordingly
repo_owner <- 'Harishneelam'
repo_name <- 'G2P-Datasets-App'
branch_name <- 'main' 
commit_message <- 'Updated Meta data via app'
access_token <- "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"

GET_meta_data <- function(){
  file_path <- 'testing-main/Meta_data_Tags.csv'
  csv_content <- get_file_from_github(repo_owner, repo_name, file_path, branch_name, access_token)
  df <- read_csv(csv_content)
  return(df)
}

get_file_from_github <- function(repo_owner, repo_name, file_path, branch_name, access_token) {
  url <- paste0('https://api.github.com/repos/', repo_owner, '/', repo_name, '/contents/', file_path, '?ref=', branch_name)
  headers <- add_headers(
    Authorization = paste('token', access_token),
    `Accept` = 'application/vnd.github.v3.raw'
  )
  response <- GET(url, headers)
  if (status_code(response) == 200) {
    return(content(response, "text"))
  } else {
    stop(paste('Error fetching file:', content(response, "text")))
  }
}

# Update Meta data in Github
update_meta_data <- function(new_data,Meta_Data){
  Meta_Data <- rbind(Meta_Data,new_data)
  update_file_on_github(repo_owner, repo_name, file_path, Meta_Data, branch_name, commit_message, access_token)
}

update_file_on_github <- function(repo_owner, repo_name, file_path, content, branch_name, commit_message, access_token) {
  # Encode the content
  csv_content <- capture.output(write.csv(content, row.names = FALSE))
  csv_string <- paste(csv_content, collapse = "\n")
  file_content <- base64encode(charToRaw(csv_string))
  # Get the SHA of the file to update it
  url <- paste0('https://api.github.com/repos/', repo_owner, '/', repo_name, '/contents/', file_path, '?ref=', branch_name)
  headers <- add_headers(
    Authorization = paste('token', access_token),
    `Content-Type` = 'application/json'
  )
  response <- GET(url, headers)
  if (status_code(response) != 200) {
    stop(paste('Error fetching file:', content(response, "text")))
  }
  file_sha <- content(response)$sha
  # Create the update payload
  data <- jsonlite::toJSON(list(
    message = commit_message,
    content = file_content,
    sha = file_sha,
    branch = branch_name
  ), auto_unbox = TRUE)
  # Make the request to update the file
  url <- paste0('https://api.github.com/repos/', repo_owner, '/', repo_name, '/contents/', file_path)
  response <- PUT(url, headers, body = data)
  if (status_code(response) == 200) {
    print('Data uploaded successfully!')
  } else {
    stop(paste('Error uploading data:', content(response, "text")))
  }
}

verify_folder <- function(folder){
  file_path <- paste0('testing-main/Datasets/',folder)
  verify_data_in_github(repo_owner, repo_name, file_path, branch_name, access_token)
}

verify_data_in_github <- function(repo_owner, repo_name, file_path, branch_name, access_token){
  url <- paste0('https://api.github.com/repos/', repo_owner, '/', repo_name, '/contents/', file_path, '?ref=', branch_name)
  headers <- add_headers(
    Authorization = paste('token', access_token),
    `Accept` = 'application/vnd.github.v3.raw'
  )
  response <- GET(url, headers)
  if (status_code(response) == 200) {
    return(TRUE)
  } else if (status_code(response) == 404) {
    return(FALSE)
  } else {
    return(FALSE)
  }
}

