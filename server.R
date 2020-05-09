
library(shiny)
library(bomrang)
library(tidyverse)
library(forecast)
library(BiocManager)
options(repos = BiocManager::repositories())
library(tidyverse)
library(class)
library(cvTools)
library(ggplot2)
library(e1071)
library(pheatmap)
library(caret)
library(survival)
library(survminer)
library(Hmisc)
library(randomForestSRC)
library(glmnet)
library(GEOquery )
library(ggplot2)
library(reshape2)


# Note: please change this dir to point to the folder where your dataset is
datadir = "GSE120396/"

# Read in the files
fileNames <- list.files(datadir)
gse = c()

for(i in 1:length(fileNames)){
  temptable <- read.delim(file.path(datadir, fileNames[i]), header=TRUE)
  gse <- cbind(gse, temptable[,2])
  colnames(gse)[i] <- colnames(temptable)[2]
}

rownames(gse) = read.delim(file.path(datadir, fileNames[1]), header=TRUE)[,1]

clinical_outcome <-getGEO("GSE120396")
clinical_outcome <- clinical_outcome$GSE120396_series_matrix.txt.gz

print(clinical_outcome$characteristics_ch1.1[1:10])

rejection_status  <- clinical_outcome$characteristics_ch1.1
rejection_status <- unlist(lapply( strsplit(as.character(rejection_status), ": " ) , `[[` , 2)  )
table(rejection_status)

largevar = apply(gse, 1, var)
ind = which(largevar > quantile(largevar, 0.9))

all_stuff = c(1)

X = (t(gse[ind,]))
y = rejection_status
re_run_knn <- function(ka) {
  cvK = 5  # number of CV folds
  cv_50acc5_knn = cv_50acc5_svm = cv_50acc5_rf = c()
  cv_acc_knn = cv_acc_svm = cv_acc_rf = c()
  
  n_sim = 25 ## number of repeats
  for (i in 1:n_sim) {
    
    cvSets = cvTools::cvFolds(nrow(X), cvK)  # permute all the data, into 5 folds
    cv_acc_knn = cv_acc_svm = cv_acc_rf = c()
    
    for (j in 1:cvK) {
      test_id = cvSets$subsets[cvSets$which == j]
      X_test = X[test_id, ]
      X_train = X[-test_id, ]
      y_test = y[test_id]
      y_train = y[-test_id]
      
      ## KNN
      fit5 = class::knn(train = X_train, test = X_test, cl = y_train, k = ka)
      cv_acc_knn[j] = table(fit5, y_test) %>% diag %>% sum %>% `/`(length(y_test))
    }
    cv_50acc5_knn <- append(cv_50acc5_knn, mean(cv_acc_knn))
  } ## end for
  return(as.data.frame(cv_50acc5_knn))
}

shinyServer(function(input, output) {a
  
  test <- reactiveValues()
  
  
  fetch_data = eventReactive(input$button, {
    a <- re_run_knn(input$slider1)
    test$df <- as.data.frame(rbind(test$df, cbind.data.frame(a, k=toString(input$slider1))))
    test$recent_k = input$slider1
    return(a)
  })
  
  observeEvent(input$cum_button, {
    test$df <- data.frame()
  })
  
  output$wind_plot = renderPlot({
    fetch_data() %>% 
      ggplot() +
      geom_boxplot(aes(y=cv_50acc5_knn)) +
      labs(y = "Accuracy", 
           title=paste0("Accuracy of KNN where K = ", toString(test$recent_k))) + 
      theme(axis.text.x = element_blank(),
          plot.title = element_text(color="Black", size=14, face="bold", hjust = 0.5),
      )
  })
  
  output$cum_plot = renderPlot({
    if(length(test$df) > 0){
      ggplot(data=test$df, aes(x=k, y=cv_50acc5_knn)) +
      geom_boxplot() +
      labs(x = "K", 
           y = "Accuracy",
           title = "Accuracy of KNN at various K") + 
        theme(
              plot.title = element_text(color="Black", size=14, face="bold", hjust = 0.5)
        )
    }
      
  })
  
})