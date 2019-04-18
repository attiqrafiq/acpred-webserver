library(protr)
library(seqinr)
#library(Interpol)
#library(caret)

library(kernlab)
library(e1071)
library(data.table)

shinyServer(function(input, output, session) {
  
  # Loads the Model to memory
  filepath <- file.path("data","Model.rds")
  mod <- readRDS(filepath)
  
  
  observe({
    
    shinyjs::hide("downloadData") # Hide download button before input submission
    if(input$submitbutton>0)
      shinyjs::show("downloadData") # Show download button after input submission
  })
  
  observe({
    FASTADATA <- ''
    fastaexample <- '>ACP1
GLWSKIKEVGKEAAKAAAKAAGKAALGAVSEAV
>ACP2
GLFDIIKKIAESI
>ACP3
GLLDIVKKVVGAFGSL
>non-ACP1
MTISLIWGIAMVVCCCIWVIFDRRRRKAGEPPL
>non-ACP2
MFATPLRQPTNASGARPAVSMDGQETPFQYEITD
>non-ACP3
LLWRKVAGATVGPGPVPA
'
    
    if(input$addlink>0) {
      isolate({
        FASTADATA <- fastaexample
        updateTextInput(session, inputId = "Sequence", value = FASTADATA)
      })
    }
  })
  
  datasetInput <- reactive({
    
    inFile <- input$file1 
    inTextbox <- input$Sequence

    if (is.null(inTextbox)) {
      return("Please insert/upload sequence in FASTA format")
    } else {
      if (is.null(inFile)) {
        # Read data from text box
        x <- inTextbox
        write.fasta(sequence = x, names = names(x),
                    nbchar = 80, file.out = "text.fasta")
        x <- readFASTA("text.fasta")
        
        # Feature extraction for Testing set
        #test <- read.fasta(x, seqtype="AA", as.string = TRUE)###read data
        test <- x
        test <- test[(sapply(test, protcheck))]###check special symbol
        AACtest <- t(sapply(test, extractAAC))
        col = 20+ 2*4
        APAACtest  <- matrix(nrow = length(test), ncol = col)
        for (i in 1:length(test)){
          APAACtest[i,] = extractAPAAC(test[[i]][1],lambda = 4, w = 0.01, customprops = NULL)
        }
        
        
        Dtest = data.frame(AACtest,APAACtest)
        
        # Predicting unknown sequences
        results <- data.frame(Prediction= predict(mod,Dtest), round(predict(mod,Dtest,type="prob"),3))
        
        print(results)
      } 
      else {  
        # Read data from uploaded file
        x <- readFASTA(inFile$datapath)
        
        # Feature extraction for Testing set
        #test <- read.fasta(x, seqtype="AA", as.string = TRUE)###read data
        test <- x
        test <- test[(sapply(test, protcheck))]###check special symbol
        AACtest <- t(sapply(test, extractAAC))
        col = 20+ 2*4
        APAACtest  <- matrix(nrow = length(test), ncol = col)
        for (i in 1:length(test)){
          APAACtest[i,] = extractAPAAC(test[[i]][1],lambda = 4, w = 0.01, customprops = NULL)
        }
        
        
        Dtest = data.frame(AACtest,APAACtest)
        
        # Predicting unknown sequences
        results <- data.frame(Prediction= predict(mod,Dtest), round(predict(mod,Dtest,type="prob"),3))
        
        print(results)
      }
    }
  })
  
  output$contents <- renderPrint({
    if (input$submitbutton>0) { 
      isolate(datasetInput()) 
    } else {
      return("Server is ready for prediction.")
    }
  })
  
  output$downloadData <- downloadHandler(
    filename = function() { paste('predicted_results', '.csv', sep='') },
    content = function(file) {
      write.csv(datasetInput(), file, row.names=FALSE)
    })
  
})
