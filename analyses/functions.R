## Stores all functions used to complete analysis that are not in the SEraster package.

calculateDensity <- function(matrix.array) {
  sum(matrix.array != 0)/(dim(matrix.array)[1] * dim(matrix.array)[2])
}

rotateAroundCenter <- function(pos_orig, angle_deg) {
  ## rotate around midpoint in both x and y axes
  pos_rotated <- rearrr::rotate_2d(data = data.frame(pos_orig), degrees = angle_deg, x_col = "x", y_col = "y", origin_fn = rearrr::midrange, overwrite = TRUE)
  
  out <- as.matrix(pos_rotated[,c("x_rotated", "y_rotated")])
  colnames(out) <- c("x", "y")
  rownames(out) <- rownames(pos_orig)
  
  ## output (class = matrix array)
  return(out)
}

## input is assumed to be a data frame with logical labels for whether it is statistically significant or not
calculatePerformanceMetrics <- function(input) {
  TP = 0
  FP = 0
  TN = 0
  FN = 0
  for (i in input$gene) {
    result.pred <- input$pred[input$gene == (i)]
    result.obs <- input$obs[input$gene == (i)]
    
    if (result.pred == TRUE && result.obs == result.pred) {
      TP <- TP + 1
    } else if (result.pred == TRUE && result.obs != result.pred) {
      FP <- FP + 1
    } else if (result.pred == FALSE && result.obs == result.pred) {
      TN <- TN + 1
    } else if (result.pred == FALSE && result.obs != result.pred) {
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

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}