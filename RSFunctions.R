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
library(tm)

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
      API_Article(Title,Title)
    }
  } else{
    response <- curl_fetch_memory(DOI, handle = handle_setheaders(new_handle(), accept = "application/x-bibtex"))
    bb <- rawToChar(response$content)
    lines <- unlist(strsplit(bb, "\n"))
    author_line <- grep("author\\s*=\\s*\\{", lines, value = TRUE)
    authors <- sub(".*?\\{(.*?)\\}.*", "\\1", author_line)
    author <- gsub(",","",strsplit(authors, " and ")[[1]][1])
    API_Article(Title,paste(Title,author))
  }
}

# Proceed if any error with API occurs
API_Article <- function(Title,SearchTerm){
  return(list("<span style='color: grey;'>Please proceed with your submission.</span>", TRUE))
  tryCatch({
    return(list("<span style='color: grey;'>Please proceed with your submission.</span>", TRUE))
    #find_Article(Title, SearchTerm)
  }, error = function(e) {
    return(list("<span style='color: grey;'>Please proceed with your submission.</span>", TRUE))
  })
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
  print(response)
  result <- fromJSON(httr::content(response, "text", encoding = "UTF-8"))
  if(result$total==0){
    return(list("<span style='color: red;'>The Article couldn't be found!</span>",FALSE)) 
  } else{
    result <- result$data[[1]]
    similarity <- 1 - stringdist(Title, result$title[1], method = "jw")
    if(similarity > 0.8){
      # Get the publication link
      if(is.null(result$externalIds$DOI[1])){
        pub_link <- 'NA'
      } else{
        doi_url <- paste0("https://doi.org/",result$externalIds$DOI[1])
        pub_link <- GET(doi_url)$url
      }
      # Get the Abstract
      abstract <- ifelse(is.null(result$abstract[1]),'',sub("^[\\s:~*]+", "", result$abstract[1]))
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
  response <- str_trim(httr::content(response)$choices[[1]]$message$content)
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
  question <- "What species is the given text is about, Give its scientific name and general common name separeted by hypen with their first letters being capital and no other text:"
  askAboutText(Text, question, chatGPTApiKey)
}

# Generate Tags - chat gpt
add_tags <- function(Text){
  chatGPTApiKey <<- "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"
  question <- "Return 5 to 8 tags for the given abstract text and seperate them by semi colons, if there is no text return NA:"
  askAboutText(Text, question, chatGPTApiKey)
}

# Generate Tags - Text mining
generate_tags <- function(Text){
  stop_words <- stopwords("en")
  capitalized_stop_words <- sapply(stop_words, function(word) {paste0(toupper(substr(word, 1, 1)), substr(word, 2, nchar(word)))})
  stop_words <- c(stop_words,capitalized_stop_words,"will","within","can","show","instead","several","across")
  Text <- gsub("-","",Text)
  Text <- gsub("\\[.*?\\]|\\(.*?\\)", "", Text)
  Text <- gsub("[^a-zA-Z]", " ", Text)
  words <- unlist(strsplit(Text, "\\s+"))
  words <- words[!grepl("ly$", words)]
  words <- words[!grepl("ing$", words)]
  words <- words[!grepl("ed$", words)]
  words <- words[nchar(words) > 1]
  words <- ifelse(words %in% stop_words, "-", words)
  cleaned_text <- paste(words, collapse = " ")
  sentence <- gsub("\\s+", " ", cleaned_text)
  sentence <- gsub("[,.-;:()/]", "-", sentence)
  words <- trimws(strsplit(sentence, "-")[[1]])
  words <- words[words != ""]
  words <- words[sapply(words, function(word) count_uppercase(word) >= 3)]
  words <- sapply(words, function(word) {paste0(toupper(substr(word, 1, 1)), substr(word, 2, nchar(word)))})
  words <- unique(words)
  print(words)
  Tags <- ifelse(length(words)>8,paste(sample(words,8), collapse = "; "),paste(words, collapse = "; "))
  print(Tags)
  return(Tags)
}

count_uppercase <- function(word) {
  # Get Important one word words
  if(sum(str_count(word, " "))==0){
    matches <- gregexpr("[A-Z]", word)
    num_uppercase <- length(unlist(regmatches(word, matches)))
    return(num_uppercase)
  }
  else{
    # Get rid of lengthy words
    if(sum(str_count(word, " "))>=3){
      return(0)
    } else{
      return(100)
    }
  }
}

# Create Meta Data
create_meta_data <- function(DOI,Pub_link,Title,Abstract,Data_Link,Spec,hasPheno,nSamples,nMarkers,Code,Meta_Data){
  # If Chat GPT API doesn't work, generating something similar
  
  # Species
  result1000 <- try({
    if(Spec==""){
      Species <- strsplit(Get_Species(paste(Title,Abstract))," - ")[[1]]
    } else{
      Species <- strsplit(Get_Species(paste(Spec))," - ")[[1]]
    }
  }, silent = TRUE)
  if(inherits(result1000, "try-error")) {
    Species <- c("",paste0(toupper(substring(Spec, 1, 1)), substring(Spec, 2)))
  }
  
  # Tags
  result2000 <- try({
    # Chat GPT
    Tags <- add_tags(paste(Title,Abstract))
  }, silent = TRUE)
  if(inherits(result2000, "try-error") || length(Tags) == 0) {
    # Basic Text Mining - not perfect but something
    Tags <- generate_tags(paste(Title,Abstract))
  }
  
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
  sentence <- gsub("[^a-zA-Z\\s]", " ", tools::toTitleCase(Title))
  words = strsplit(sentence, " ")[[1]]
  unique_txt <- paste(sort(tail(words[order(nchar(words))],3)), collapse = " ")
  Unique_name <- gsub("[-',:() ]", "", paste(unique_txt, Species_Common_Name))
  Folder_name <- gsub(" ","",paste(sprintf("%05d",nrow(Meta_Data)+1),"_",Unique_name))
  return(Folder_name)
}

