datasets <- list(
  "ui" = page_fillable(
    
    layout_columns(
      width = c(6,6),
      style = "height:100vh;",
      
      card(
        fluidRow(
          column(12, DTOutput("metadataTable")) 
        ),
        style = "height=100%;"
      ),
      
      card(
           htmlOutput("datasetSummary")
    ))
  ),
  
  "server" = function(input, output, session){
    
    output$metadataTable <- renderDT(datatable(
      assembleMetadata(metadata),
      selection = "single",
      filter = "top",
      rownames = F,
      options = list(
        pageLength = -1,
        deferRender = T,
        scrollY = TRUE,
        scrollX = TRUE,
        dom = 'ft',
        autoWidth = TRUE,
        responsive = TRUE
      ),
      
      # Search functionality
      callback = JS(
        "$(document).ready(function() {",
        "  $('#example').DataTable();",
        "  $.fn.dataTable.ext.search.push(",
        "    function(settings, searchData, index, rowData, counter) {",
        "      var match = false;",
        "      var searchTerm = settings.oPreviousSearch.sSearch.toLowerCase();",
        "      searchData.forEach(function(item, index) {",
        "        if (item.toLowerCase().startsWith(searchTerm)) {",
        "          match = true;",
        "        }",
        "      });",
        "      return match;",
        "    }",
        "  );",
        "});"
      ),
      escape = F
    ))
    
    # Default selection
    observe({
      DT::selectRows(dataTableProxy("metadataTable"), 1)
    })
    
    # User action
    observeEvent(input$metadataTable_rows_selected, {
      
      selected_row <- input$metadataTable_rows_selected
      selectedData <- metadata[selected_row,]
      
      output$datasetSummary <- renderUI({
        HTML(paste0("<h3 style='margin-bottom: 10px;'>", selectedData$Title, "</h3>",
                    
                    "<p style='margin-bottom: 10px;'><strong>DOI:</strong> <a href='",
                    selectedData$DOI, "' target='_blank'>", selectedData$DOI, "</a></p>",
                    
                    "<p style='margin-bottom: 20px;'><strong>Scientific name, Common name:</strong> ", 
                    selectedData$Scientific_Name, ", ", selectedData$Common_Name, "</p>",
                    
                    accordion(
                      id = "Abstract",
                      open = FALSE,
                      accordion_panel(
                        title = HTML("<strong>Abstract</strong>"),
                        HTML(paste0("<p>", selectedData$Abstract, "</p>"))
        )),
        
        accordion(
          id = "code",
          open = FALSE,
          accordion_panel(
            title = HTML("<strong>Coad to load, curate, and analyze</strong>"),
            
            HTML(paste0("<p>To download the dataset, go to <a href='", selectedData$Sharing_Link, "' target='_blank'>",
                        selectedData$Sharing_Link, "</a> and follow download instructions.</p>")),
            HTML("<p>To load, curate, and analyze the dataset, create an R session and set the working directory to the downloaded data folder.</p>"),
            HTML("<br>"),
            
            HTML(paste0("<p>", "To load data, run this code:", "</a></p>")),
            verbatimTextOutput("loaderCode", placeholder = TRUE),
            uiOutput("copyButtonLoader"),
            HTML("<br>"),
            
            HTML(paste0("<p>", "To curate data, run this code:", "</a></p>")),
            verbatimTextOutput("curationCode", placeholder = TRUE),
            uiOutput("copyButtonCuration"),
            HTML("<br>"),
            
            HTML(paste0("<p>", "To analyze data, run code below:", "</a></p>")),
            navset_pill(
              id = "methodTabs",
              !!!lapply(method_names, function(method) {
                nav_panel(
                  title = method,
                  verbatimTextOutput(paste0("analysisCode_", method), placeholder = TRUE),
                  uiOutput(paste0("copyButtonAnalysis_", method)),
                  tags$br()
                )
              })
            )
          )),
        
        accordion(
          id = "Citation",
          open = FALSE,
          accordion_panel(
            title = HTML(paste0("<p style='font-style: italic'>", 
                                "<strong>Citation information</strong>")),
            verbatimTextOutput("Citation", placeholder = TRUE),
            uiOutput("copyButtonCitation"),
          )
        )
        ))
      })
      
      
      loaderCode <- assembleLoaderCode(selected_row)
      output$copyButtonLoader <- renderUI({
        rclipButton(
          inputId = "copy", 
          label = "",
          clipText = loaderCode,
          icon = icon("clipboard"),
          style = "padding: 5px 10px"
        )
      })
      output$loaderCode <- renderText(loaderCode)
      
      
      curationCode <- assembleCurationCode(selected_row)
      output$copyButtonCuration <- renderUI({
        rclipButton(
          inputId = "copy", 
          label = "",
          clipText = curationCode,
          icon = icon("clipboard"),
          style = "padding: 5px 10px"
        )
      })
      output$curationCode <- renderText(curationCode)
      
      
      # analysisCode <- assembleAnalysesCode(selected_row)
      # output$copyButtonAnalysis <- renderUI({
      #   rclipButton(
      #     inputId = "copy", 
      #     label = "",
      #     clipText = analysisCode,
      #     icon = icon("clipboard"),
      #     style = "padding: 5px 10px"
      #   )
      # })
      # output$analysisCode <- renderText(analysisCode)
      
      
      # Define the directory containing method files
      AnalysesCodeDir <- "GPDatasets-main/Analyses"
      
      # Dynamically generate method names from file names
      method_files <- list.files(AnalysesCodeDir, pattern = "\\.R$", full.names = TRUE)
      method_names <- gsub("\\.R$", "", basename(method_files))
      
      # Read and prepare code for each method
      method_codes <- lapply(method_files, readLines)
      names(method_codes) <- method_names
      
      # Render code and copy buttons for each method dynamically
      lapply(seq_along(method_names), function(i) {
        method <- method_names[i]
        
        # Render the analysis code
        output[[paste0("analysisCode_", method)]] <- renderText({
          paste(method_codes[[i]], collapse = "\n")
          })
        
        # Render the copy button for the analysis code
        output[[paste0("copyButtonAnalysis_", method)]] <- renderUI({
          rclipButton(
            inputId = paste0("copy_", method),
            label = "",
            clipText = paste(method_codes[[i]], collapse = "\n"),
            icon = icon("clipboard"),
            style = "padding: 5px 10px"
          )
        })
      })
      
      # Assemble citation
      Citation <- assembleDatasetCitation(selected_row)
      output$copyButtonCitation <- renderUI({
        rclipButton(
          inputId = "copy", 
          label = "",
          clipText = Citation,
          icon = icon("clipboard"),
          style = "padding: 5px 10px"
        )
      })
      output$Citation <- renderText(Citation)
      
    })
    
  }
  
)

