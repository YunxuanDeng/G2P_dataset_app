text1 <- "The G2P-Datasets app is a platform for accessing public genome-to-phenome datasets. The app provides searchable metadata on >100 plant and animal datasets paired with code to load, curate, and analyze data. Users can contribute new datasets to expand G2P-Datasets as a resource for benchmarking, data discovery, and education. G2P-Datasets is documented at (placeholder for publication DOI), and the open-source app code and metadata are available at (a placeholder for final GitHub repo)."
text2 <- "G2P-Datasets is a web app that requires only a browser. To browse genome-to-phenome datasets and access code to load, curate, and analyze data, go to the \"Browse datasets\" module. To contribute another genome-to-phenome dataset and/or analysis to G2P-Datasets, go to the \"Add dataset\" module."
text3 <- "Please follow the instructions carefully in the following document to submit a new dataset."

about <- list(
  "ui" = page_fillable(htmlOutput("about_info")),
  
  "server" = function(input, output, session){
    output$about_info = renderUI({
      HTML(paste0("<h5><strong>", "What is G2P-Datasets?", "</strong></h5>",
                  "<p>", text1, "</p>",
                  "<br>",
                  "<h5><strong>", "How do I use G2P-Datasets?", "</strong></h5>",
                  "<p>", text2, "</p>",
                  "<br>",
                  "<h5><strong>", "How do I submit a dataset?", "</strong></h5>",
                  "<p>", text3, "</p>",
                  "<a href='https://drive.google.com/file/d/1NeJ3CxIhUVJZEGCtVmuXZq7j2lYvB5jw/view?usp=sharing' target='_blank'>",
                  "📄 Instructions to submit a dataset","</a>",
                  "<br>",
                  "<br>",
                  "<p>", "To report erroneous data or bugs in the app, please contact <a href=\"mtwason@ucdavis.edu\"> mtwaston@ucdavis.edu </a>.", "</p>",
                  "<br>",
                  "<img src='https://storage.googleapis.com/kaggle-organizations/4421/thumbnail.png?t=2025-01-17-04-03-19' width='200' />"))
      })
  }
)