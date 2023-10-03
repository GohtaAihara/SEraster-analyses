## Stores all functions used to complete analysis that are not in the SEraster package.

calculateDensity <- function(matrix.array) {
  sum(matrix.array != 0)/(dim(matrix.array)[1] * dim(matrix.array)[2])
}

## input is assumed to be a data frame with logical labels for whether it is statistically significant or not
calculatePerformanceMetrics <- function(input) {
  TP = 0
  FP = 0
  TN = 0
  FN = 0
  for (i in input$gene) {
    result.px <- input$pixel[input$gene == (i)]
    result.sc <- input$singlecell[input$gene == (i)]
    
    if (result.px == TRUE && result.sc == result.px) {
      TP <- TP + 1
    } else if (result.px == TRUE && result.sc != result.px) {
      FP <- FP + 1
    } else if (result.px == FALSE && result.sc == result.px) {
      TN <- TN + 1
    } else if (result.px == FALSE && result.sc != result.px) {
      FN <- FN + 1
    }
    # matrix(c(TP,FN,FP,TN), nrow = 2, ncol = 2)
    TPR <- TP/(TP + FN)
    FPR <- FP/(FP + TN)
    specificity <- 1 - FPR
    PPV <- TP/(TP + FP)
    F1 <- 2 * (TPR * PPV) / (TPR + PPV)
    ACC <- (TP + TN) / (TP + FP + TN + FN)
  }
  return(data.frame(TP=TP,FP=FP,TN=TN,FN=FN,TPR=TPR,FPR=FPR,specificity = specificity, PPV=PPV,F1=F1,ACC=ACC))
}
