modulePanel <- function(title, value){
  tabPanel(
    title = title,
    value = value,
    eval(parse(text = value))$ui
  )
}

reloadWarning <- tags$head(tags$script(HTML("
    // Enable navigation prompt
    window.onbeforeunload = function() {
        return 'Your changes will be lost!';
    };
")))

unescape_html <- function(str){
  xml2::xml_text(xml2::read_html(paste0("<x>", str, "</x>")))
}

readURL <- function(url){
  req <- GET(url)
  stop_for_status(req)
  return(content(req))
}

askChatGPT <- function(prompt, chatGPTApiKey){
  
  response <- httr::POST(
    url = "https://api.openai.com/v1/chat/completions", 
    add_headers(Authorization = paste("Bearer", chatGPTApiKey)),
    content_type("application/json"),
    encode = "json",
    body = list(
      model = "gpt-3.5-turbo",
      messages = list(
        list(role = "user", content = prompt)
      ),
      temperature = 0.7,
      max_tokens = 512
    )
  )

  response <- str_trim(content(response)$choices[[1]]$message$content)
  return(response)
  
}

getSummary <- function(doi, chatGPTApiKey){
  response <- paste0("http://doi.org/", doi) %>%
    readURL() %>%
    unescape_html() %>%
    substr(1, 1000) %>%
    paste0("Summarise this publication in one sentence: ", .) %>%
    askChatGPT(., chatGPTApiKey)
  return(response)
}

