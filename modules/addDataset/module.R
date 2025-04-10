addDataset <- list(
  
  "ui" = fluidPage(
    useShinyjs(),
    h1("Add dataset"),
    p("The GPDatasets repository is designed so that users can add their own public datasets. To add a dataset/database, first fill in the fields below, download the generated dataset folder, then push the dataset folder to", a("GPDatasets GitHub", href = "https://github.com/QuantGen/GPDatasets")),
    HTML("<strong>Important note</strong><span style='color:red;'>*</span>: The data should atleast contain genotypic data(SNPs) for a species."),
    tags$head(
      tags$style(HTML("
    .input-group {
        display: flex;
        flex-direction: column; 
      }
      .input-group .input-row {
        display: flex;
        align-items: center; 
      }
      .input-group .form-control {
        flex: 2;
      }
      .input-group .btn {
        margin-left: 10px;
        margin-top: 10px;
      }
      .grey-text {
        color: grey;
        font-size: 0.9em; 
        margin-top: 0px; 
      }
      
    }
    "))
    ),
    br(),
    br(),
    fluidRow(
      column(6,
             div(class = "input-group", 
                 div(class = "input-row",
                     textInput("addDOI", HTML("<strong>Data DOI</strong><span style='color:red;'>*</span><span style='font-weight:normal; color:grey; font-size:12px;'> ( Example: https://doi.org/XYZ )</span>"), width = "80%", value = "https://doi.org/10.3389/fpls.2021.803736"),
                     actionButton("validateDOI", "Validate DOI", width = "20%")),
                 tags$span("Note: If the data is part of an article and included within it, 
                       or if the external source hosting the data lacks a DOI, include the article's DOI 
                       and provide a direct sharing link to the data.", class = "grey-text")
             )),
      column(6,
             div(class = "input-group",
                 div(class = "input-row",
                     textInput("addLink", HTML("<strong>Data sharing link</strong><span style='color:red;'>*</span>"), width = "80%", value = ""),
                     actionButton("validateLink", "Validate Link", width = "20%")
                 )))),
    fluidRow(
      column(6,uiOutput("result1")),
      column(6,htmlOutput("result2"))),
    
    fluidRow(column(12,DTOutput("resultTable"))),
    br(),
    
    fluidRow(
      column(8,
             div(class = "input-group",
                 div(class = "input-row",
                     textInput("addTitle", HTML("<strong>Publication Title</strong><span style='color:red;'>*</span><span style='font-weight:normal; color:grey; font-size:12px;'> ( Enter the title of the article as published. )</span>"),width = "80%",value = "CottonGVD: A Comprehensive Genomic Variation Database for Cultivated Cottons"),
                     actionButton("findArticle", "Find Article",width = "20%"))
             )),
      column(4,
             div(class = "input-group",
                 div(class = "input-row",
                     textInput("addSpec", HTML("<strong>Species</strong>"),width = "100%",value = "")
                 ))),
      
    ),
    fluidRow(column(8,uiOutput("result3"))),
    fluidRow(column(8,uiOutput("modifyArticle"))),
    fluidRow(
      column(5,
             tags$div(
               HTML("<strong>Code to load data</strong><span style='color:grey; font-size:12px;'> ( Refer to existing datasets )</span>"),
               style = "margin-bottom: 5px;"
             ),
             textAreaInput("addLoaderCode",NULL, height = "110px", width = "100%",
                           value = "# Recommended
# Import necessary packages 

# Get geno
geno <- ''"),
             checkboxInput("hasPheno", "Data includes Phenotypic information", value = FALSE)),
      
      fluidRow(column(3, 
                      numericInput("nSamples",  HTML("<strong>Number of samples</strong><span style='color:red;'>*</span>"), value = 0, min = 0, step = 1),
                      numericInput("nMarkers",  HTML("<strong>Number of markers(SNPs)</strong><span style='color:red;'>*</span>"), value = 0, min = 0, step = 1))),
    ),
    br(),
    
    fluidRow(column(3,actionButton("createData", "Preview"))),
    fluidRow(column(8,uiOutput("result4"))),
    br(),
    fluidRow(column(12,DTOutput("metaData"))),
    br(),
    fluidRow(column(8,uiOutput("proceedButton"))),
    br(),
    verbatimTextOutput("status"),
    fluidRow(column(8,uiOutput("result5"))),
    br(),
    fluidRow(column(8, uiOutput("result6"))),
    fluidRow(column(6, downloadButton("downloadDataset", "Download dataset folder",style = "margin-left: 22px;"))),
    fluidRow(column(8, uiOutput("result7"))),
    rclipboardSetup(),
    fluidRow(column(5,uiOutput("Code1"),style = "padding: 0; margin-left: 22px;"),
             column(1,uiOutput("copyButton1"),style = "padding: 0;margin-left: 0;")),
    fluidRow(column(5,uiOutput("Code2"),style = "padding: 0; margin-left: 22px;"),
             column(1,uiOutput("copyButton2"),style = "padding: 0;margin-left: 0;")),
    fluidRow(column(5,uiOutput("Code3"),style = "padding: 0; margin-left: 22px;"),
             column(1,uiOutput("copyButton3"),style = "padding: 0;margin-left: 0;")),
    fluidRow(column(5,uiOutput("Code4"),style = "padding: 0; margin-left: 22px;"),
             column(1,uiOutput("copyButton4"),style = "padding: 0;margin-left: 0;")),
    fluidRow(column(5,uiOutput("Code5"),style = "padding: 0; margin-left: 22px;"),
             column(1,uiOutput("copyButton5"),style = "padding: 0;margin-left: 0;")),
    fluidRow(column(5,uiOutput("Code6"),style = "padding: 0; margin-left: 22px;"),
             column(1,uiOutput("copyButton6"),style = "padding: 0;margin-left: 0;")),
    fluidRow(column(5,uiOutput("Code7"),style = "padding: 0; margin-left: 22px;"),
             column(1,uiOutput("copyButton7"),style = "padding: 0;margin-left: 0;")),
    fluidRow(column(8,uiOutput("result8"))),
    fluidRow(column(8,uiOutput("submitButton"))),
    br()
  ),
  
  "server" = function(input, output, session){
    ## Get the metadata from GitHub repo
    Meta_Data <- metadata
    # Meta_Data <- GET_meta_data() # update final repo details in RSFunction
    # Meta_Data <- readxl::read_xlsx('/Users/harishneelam/Desktop/Quantgen/GPDatasets/Meta_data.xlsx')
    
    Bresult1 <- reactiveVal(NULL)
    observeEvent(input$validateDOI, {
      # Checks if the doi is null, already exists in the data and whether it is a valid url.
      showNotification("Validating DOI...", duration = NULL, id = "processNotificationdoi")
      result1 <- check_url(input$addDOI,Meta_Data) 
      Bresult1(result1)
      removeNotification("processNotificationdoi")
      output$result1 <- renderText({
        HTML(paste(result1[[1]]))
      })
      if(result1[[2]][1]==TRUE){
        # Fill the Data sharing link automatically as both are connected most often
        updateTextInput(session, "addLink", value = result1[[3]][1])
        
        # Display the existing data if exists already
        if(result1[[4]][1]==1){
          output$resultTable <- renderDT({
            datatable(result1[[5]], options = list(
              pageLength = 10,
              searching = FALSE,
              lengthChange = FALSE,
              paging = FALSE,
              info = FALSE
            ),selection = "none")
          })
        }
      } 
      
    })
    
    
    Bresult2 <- reactiveVal(NULL)
    observeEvent(input$validateLink, {
      # Inputs the sharing from doi and checks if the link is valid.
      showNotification("Validating Link...", duration = NULL, id = "processNotificationlink")
      result2 <- check_link(input$addLink)
      Bresult2(result2)
      removeNotification("processNotificationlink")
      output$result2 <- renderText({
        HTML(result2[[1]][1])
      })
    })
    
    Bresult3 <- reactiveVal(NULL)
    observeEvent(input$findArticle, {
      # finds the article based on Title and author information obtained from doi(if given)
      showNotification("Finding article...", duration = NULL, id = "processNotification")
      tryCatch({
        result3 <- check_Title(input$addTitle,input$addDOI,check_url(input$addDOI,Meta_Data)[[2]][1])
        Bresult3(result3)
        removeNotification("processNotification")
        output$result3 <- renderUI({
          HTML(paste(result3[[1]]))
        })
      }, error = function(e) {
        removeNotification("processNotification")
        output$result3 <- renderUI({
          HTML("Couldn't find the article, please try refreshing the page!")
        })
      })
      # In some cases, the data doesn't belong to an article, maybe of an database, package or external source, in such cases, we let users override the request and proceed
      if(result3[[2]][1]==FALSE){
        output$modifyArticle <- renderUI({
          checkboxInput("modifyArticle", "Not an Article?",value = FALSE)
        })
      }
      
    })
    
    
    # Create meta data
    Bresult4 <- reactiveVal(NULL)
    rv <- reactiveVal(NULL)
    statusMessage <- reactiveVal(NULL)
    observeEvent(input$createData, {
      
      SBresult1 <- Bresult1()
      SBresult2 <- Bresult2()
      SBresult3 <- Bresult3()
      tryCatch({
        # Only if the details are valid
        if (as.logical(SBresult1[2]) && as.logical(SBresult2[2]) && (as.logical(SBresult3[2]) || input$modifyArticle) ) {
          withProgress(message = 'Processing...', {
            incProgress(1/2, detail = "Please wait")
            # Creates metadata
            result41 <- create_meta_data(input$addDOI,SBresult3[3],input$addTitle,SBresult3[4],
                                         input$addLink,input$addSpec,input$hasPheno,input$nSamples,input$nMarkers,
                                         input$addLoaderCode,Meta_Data)
            metaData <- cbind(result41[[1]], result41[[2]])
            colnames(metaData) <- c("Input","Value")
            rv(metaData)
          })
          
          # Checking if data exists already in case of previous message "Data exists already!" 
          res <- metaData[,"Value"]
          names(res) <- colnames(Meta_Data)
          # comparing doi, species names, n geno, n markers
          res <- res[c(2,7,13,14)]
          res[c(1,2)] <- as.character(res[c(1,2)])
          res[c(3,4)] <- as.double(res[c(3,4)])
          res <- as.data.frame(res)[1,]
          mtd <- as.data.frame(Meta_Data[,c(2,7,13,14)])
          
          if(SBresult1[4]==1 && (tail(duplicated(rbind(mtd,res)),1)>0)){
            output$result4 <- renderUI({
              tags$div(
                tags$div(style = "color: red;", "Data exists!"),
                tags$div(style = "color: black;", "Data matches with the existing data. Please try adding different data.")
              )
            })
          } else{
            output$result4 <- renderUI({
              tags$div(
                tags$div(style = "color: green;", "Data processed!"),
                tags$div(style = "color: black;", "Please verify the details below. If necessary, modify them by double-clicking the fields and fill in any missing information to make the data more comprehensive.")
              )
            })
            output$metaData <- renderDT({
              # Showing metadata
              datatable(metaData, options = list(
                pageLength = 15,
                searching = FALSE,
                lengthChange = FALSE,
                paging = FALSE,
                info = FALSE,
                columnDefs = list(
                  list(width = '150px', targets = 0), 
                  list(width = '300px', targets = 1)  
                )
              ), editable = list(target = 'cell', disable = list(columns = c(0))),selection = "none")
            })
            observeEvent(input$metaData_cell_edit, {
              current_data <- rv()
              edited_value <- input$metaData_cell_edit$value
              row_index <- input$metaData_cell_edit$row 
              # Restrict editing 
              if (row_index %in% c(2, 4, 8, 10, 13, 14)) {
                showModal(modalDialog(
                  title = "Alert!",
                  "Editing this value is not allowed.",
                  easyClose = TRUE,
                  footer = NULL
                ))
                output$metaData <- renderDT({
                  # Showing metadata
                  datatable(current_data, options = list(
                    pageLength = 15,
                    searching = FALSE,
                    lengthChange = FALSE,
                    paging = FALSE,
                    info = FALSE,
                    columnDefs = list(
                      list(width = '150px', targets = 0), 
                      list(width = '300px', targets = 1)
                    )
                  ), editable = list(target = 'cell', disable = list(columns = c(0))),selection = "none")
                })
              } else {
                new_data <- current_data
                new_data[row_index, "Value"] <- edited_value
                rv(new_data)
              }
            })
            
            output$proceedButton <- renderUI({
              actionButton("proceedData", "Proceed")
            })
          } 
        }
        else {
          output$result4 <- renderUI({
            tags$div(style = "color: red;", paste("Error! Please ensure the required fields are filled."))
          })
        }
      }, error = function(e) {
        output$result4 <- renderUI({
          tags$div(style = "color: red;", paste("Error! Please ensure the required fields are valid."))
        })
      })
    })
    
    
    
    
    dataSubmitted <- reactiveVal(FALSE)
    observeEvent(input$proceedData, {
      dataSubmitted(TRUE)
      final_data <- rv()
      final_data <- final_data[,"Value"]
      names(final_data) <- colnames(Meta_Data)
      final_data[1:8] <- as.character(final_data[1:8])
      final_data[9:12] <- as.logical(final_data[9:12])
      final_data[15:16] <- as.character(final_data[15:16])
      Bresult4(final_data)
      output$result5 <- renderUI({
        tags$div(style = "color: green;", paste("Data verified!!"))
      })
    })
    
    
    # Enable/disable the download button based on data submission
    observe({
      shinyjs::toggleState("downloadDataset", dataSubmitted())
    })
    
    # Instructions to push data to GitHub Repository
    observeEvent(input$proceedData, {
      output$result6 <- renderUI({
        tagList(
          tags$div(style = "color: black; font-weight: bold;", "Instructions to add the dataset to the GitHub Repository:"),
          tags$ul(
            tags$li("Click the following button to download the dataset folder."),
            tags$li("Ensure the file is saved to your DOWNLOADS folder.")
          )
        )
      })
      output$result7 <- renderUI({
        tagList(
          tags$ul(
            tags$li("Open your terminal and paste the following commands. Modify them as needed")
          )
        )
      })
      code_text1 <- "cd"
      output$Code1 <- renderUI({tags$pre(tags$code(code_text1))})
      output$copyButton1 <- renderUI({rclipButton(inputId = "copy",label = "",clipText = code_text1,icon = icon("clipboard"))})
      code_text2 <- "git clone git@github.com:QuantGen/GPDatasets.git "
      output$Code2 <- renderUI({tags$pre(tags$code(code_text2))})
      output$copyButton2 <- renderUI({rclipButton(inputId = "copy", label = "",clipText = code_text2,icon = icon("clipboard"))})
      code_text3 <- "cd GPDatasets/Datasets"
      output$Code3 <- renderUI({tags$pre(tags$code(code_text3))})
      output$copyButton3 <- renderUI({rclipButton(inputId = "copy", label = "",clipText = code_text3,icon = icon("clipboard"))})
      code_text4 <- paste0("mv ~/Downloads/",Bresult4()[[1]])
      output$Code4 <- renderUI({tags$pre(tags$code(code_text4))})
      output$copyButton4 <- renderUI({rclipButton(inputId = "copy", label = "",clipText = code_text4,icon = icon("clipboard"))})
      code_text5 <- "git add ."
      output$Code5 <- renderUI({tags$pre(tags$code(code_text5))})
      output$copyButton5 <- renderUI({rclipButton(inputId = "copy", label = "",clipText = code_text5,icon = icon("clipboard"))})
      code_text6 <- paste0("git commit -m 'Added ",Bresult4()[[1]], "'")
      output$Code6 <- renderUI({tags$pre(tags$code(code_text6))})
      output$copyButton6 <- renderUI({rclipButton(inputId = "copy", label = "",clipText = code_text6,icon = icon("clipboard"))})
      code_text7 <- "git push git@github.com:QuantGen/GPDatasets.git"
      output$Code7 <- renderUI({tags$pre(tags$code(code_text7))})
      output$copyButton7 <- renderUI({rclipButton(inputId = "copy", label = "",clipText = code_text7,icon = icon("clipboard"))})
      
      output$result8 <- renderUI({
        url <- "https://github.com/QuantGen/GPDatasets"  
        tagList(
          tags$ul(
            tags$li(
              "Please verify if the Dataset is uploaded at ",
              tags$a(href = url, target = "_blank", "GPDatasets GitHub")
            )
          )
        )
      })
      
      # Final submit button
      output$submitButton <- renderUI({
        actionButton("submitData", "Submit")
      })
    })
    
    
    
    observeEvent(input$submitData, {
      # Final submission
      
      # Verify if the data has been pushed and then update the meta data.
      result_final <- verify_folder(Bresult4()[[1]])
      if(result_final == TRUE){
        ################### Update Meta data csv file in Github
        # update_meta_data(Bresult4(),Meta_Data)
        showModal(modalDialog(
          title = "Thank you!",
          "The data has been submitted.",
          "Please refresh the page to add a new dataset.",
          easyClose = TRUE,
          footer = NULL
        ))
      } else{
        showModal(modalDialog(
          title = "Error!",
          "The data has not been submitted.",
          "Please follow the instructions carefully and verify if the data has been pushed.",
          easyClose = TRUE,
          footer = NULL
        ))
      }
      
    })
    
    # Generate Data Folder to download
    output$downloadDataset <- downloadHandler(
      filename = function() {
        new_meta_data <- Bresult4()
        paste0(new_meta_data[1],".zip")
      },
      content = function(file) {
        new_meta_data <- Bresult4()
        
        # Create a temporary directory for the files
        tempdir <- tempdir()
        temp_folder_path <- file.path(tempdir, new_meta_data[[1]])
        dir.create(temp_folder_path, showWarnings = FALSE, recursive = TRUE)
        original_wd <- getwd()
        on.exit(setwd(original_wd))
        setwd(temp_folder_path)
        
        # read_data_code.R generator
        writeLines(strsplit(input$addLoaderCode, "\n")[[1]][-1], "read_data_code.R")
        # create output folder
        output_folder <- "output"
        if (!dir.exists(output_folder)){
          dir.create(output_folder)
        }
        setwd(output_folder)
        # JSON generator
        writeLines(jsonlite::toJSON(new_meta_data, pretty = TRUE, auto_unbox = TRUE), "meta_data.json")
        # .md file generator
        md_content <- paste0(
          "# ", new_meta_data[4], "\n\n",
          "This publication discusses: ", new_meta_data[5], "\n\n",
          "It contains ", new_meta_data[13], " genotypes and ", new_meta_data[14], " markers.\n\n",
          "Title: ", new_meta_data[4], "\n\n",
          "Scientific name: ", new_meta_data[6], "\n\n",
          "Common name: ", new_meta_data[7], "\n\n",
          "DOI: ", new_meta_data[2], "\n\n"
        )
        writeLines(md_content, "meta_data.md")
        # .bib file generator
        url <- as.character(new_meta_data[2])
        try({
          file_name <- "citation.bib"
          curl_download(url, destfile = file_name, handle = handle_setheaders(new_handle(), accept = "application/x-bibtex"))
        })
        # zip file
        zip::zip(
          zipfile = file,
          files = dir(temp_folder_path),
          root = temp_folder_path
        )
        unlink(tempdir)
        
      },
      contentType = "application/zip"
    )
    
    
    
  }
  
  
)
