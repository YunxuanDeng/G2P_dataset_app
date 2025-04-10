require(shiny)
require(httr)
require(RCurl)
require(jsonlite)
require(bslib)
require(data.table)
require(DT)
require(rclipboard)
require(tidyverse); require(magrittr)
library(stringr) 
library(shinyjs)
library(gert)
source("utils.R")

repoDir <- "https://github.com/QuantGen/GPDatasets/archive/refs/heads/main.zip"
download.file(url = repoDir, destfile = "repo.zip")
unzip(zipfile = "repo.zip")

datasetDir <- "GPDatasets-main"

metadata <<- datasetDir %>%
  list.files(full.names = TRUE, recursive = TRUE) %>%
  {.[map_lgl(., function(x){grepl("meta_data.json", x)})]} %>%
  map(., function(x){jsonlite::fromJSON(x)}) %>%
  map(., as.data.frame) %>%
  map(., function(x){x %>% mutate(across(everything(), as.character))}) %>%
  bind_rows()

# save(metadata, file = "environment.RData")
# load("environment.RData")

sapply(list.files(path = "modules", recursive = TRUE, pattern = "^.*\\.R$", full.names = TRUE), source)
