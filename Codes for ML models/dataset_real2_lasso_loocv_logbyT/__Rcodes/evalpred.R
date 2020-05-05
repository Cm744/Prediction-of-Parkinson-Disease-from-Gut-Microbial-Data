
evaluate_pred <- function (probs_pred, y, caseid = names (y), showplot = FALSE, method = "Prediction")
{
  
  probs_attrue_bplr <- function (probs_pred, y)
  {
    tp <- rep(0, nrow(probs_pred))
    for(i in 1:nrow(probs_pred)) tp[i] <- probs_pred[i,y[i]]
    
    tp
  }
  
  ##############################################################################
  
  
  eval_tab_pred <- function (table_eval, showplot = TRUE, method = "Prediction", ...)
  {
    if (is.character (table_eval)) 
    {
      table_eval <- as.matrix (read.table (table_eval))
    }
    
    C <- ncol (table_eval) - 3
    colnames (table_eval) <- c("Case ID", "True Label", paste ("Pred. Prob", 1:C), "Wrong?")
    
    probs_pred <- table_eval [, 2+(1:C)]
    y <- table_eval[,2]
    probs_at_truelabels <- probs_attrue_bplr (probs_pred, y)
    which.wrong <- which (table_eval[,C+3] == 1)
    n <- nrow (table_eval)
    
    amlp <- - mean (log (probs_at_truelabels))
    no_errors <- sum (table_eval[, C+3])
    er <- no_errors/n
    
    pred_roc <- NA
    pred_auc <- NA
    if (C==2) {
      pred_roc <- roc(y, probs_pred[,2], quiet = TRUE)
      pred_auc <- auc(pred_roc)
    }
    
    yl <- y; if (C == 2) yl[y==2] <- 3
    
    plotargs <- list (...)
    if (is.null (plotargs$ylab)) 
      plotargs$ylab <- "Predictive Probability at True Label"
    if (is.null (plotargs$xlab)) plotargs$xlab <- "Case Index"
    if (is.null (plotargs$ylim)) plotargs$ylim <- c(0,1)
    if (is.null (plotargs$pch)) plotargs$pch <- yl   
    if (is.null (plotargs$col)) plotargs$col <- 1+table_eval[, C+3]   
    if (showplot)
    {
      plotargs$x <- probs_at_truelabels
      do.call (plot, plotargs)   
      
      if (C == 2) abline (h = 0.5)

      title (main = sprintf ("AMLP = %5.3f, Error Rate = %4.2f%% (%d/%d), AUC=%5.3f", 
                             amlp, er*100, no_errors, n, pred_auc), 
             cex = 0.8, line = 0.5) 
      
      # if (no_errors > 0) {   
      #	    text (which.wrong, probs_at_truelabels[which.wrong], labels = which.wrong,
      #             srt = 90, adj = - 0.4, cex = 0.9, col = "red")
      
      # }
    }
    pred0 <- table(y)/length(y)
    er0 <- min(pred0)
    amlp0 <- -sum(pred0*log(pred0))
    list (probs_at_truelabels = probs_at_truelabels, table_eval = table_eval,  
          amlp = amlp, R2_amlp = (amlp0-amlp)/amlp0, er = er, R2_er = (er0-er)/er0, 
          roc = pred_auc, auc = pred_auc, which.wrong = which.wrong )
  }
  
  ##############################################################################
  
  if (is.null (caseid)) caseid <- 1:length (y)
  
  C <- ncol (probs_pred)
  values_pred <- apply (probs_pred, 1, which.max)
  
  table_eval <- data.frame (caseid, y, probs_pred, 1 * (values_pred != y)) 
  colnames (table_eval) <- c("Case ID", "True Label", paste ("Pred. Prob", 1:C), "Wrong?")
  
  eval_tab_pred (table_eval, showplot = showplot, method = method)
}


