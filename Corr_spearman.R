Correlation_Pairwise_2 <- function (x,y )
{
  Corr_Coefficient_Results <- data.frame()
  Corr_P_Value_Results <- data.frame()
  
  CLASS <- colnames(x)
  for (i in seq_along(CLASS))
  {
    CORRELATION <- corr.test(x = data.frame(x[,i]), y = y, use = "pairwise", method = "spearman", adjust = "fdr", ci = FALSE)
    CORRELATION_COEFFICIENT <- CORRELATION[[1]]
    CORRELATION_P_VALUE <- CORRELATION[[4]]
    Param_Col <- colnames(CORRELATION_COEFFICIENT)
    print (i)
    for (j in seq_along(Param_Col))
    {
      #if(CORRELATION_P_VALUE[1,j] <= 1)
      #{
      newline_coefficient_j <- data.frame(t(c(paste(CLASS[[i]]),colnames(CORRELATION_COEFFICIENT)[j], CORRELATION_COEFFICIENT[1,j])))
      newline_p_value_j <- data.frame(t(c(paste(CLASS[[i]]), colnames(CORRELATION_P_VALUE)[j], CORRELATION_P_VALUE[1,j])))
      Corr_Coefficient_Results <- rbind(Corr_Coefficient_Results, newline_coefficient_j)
      Corr_P_Value_Results <- rbind(Corr_P_Value_Results, newline_p_value_j)
      #}
    }
  }
  colnames(Corr_Coefficient_Results) <- c("Taxa", "Path", "Correlation_Coefficient")
  colnames(Corr_P_Value_Results) <- c("Taxa", "Path", "Correlation_P_Value")
  write.table(Corr_Coefficient_Results, file = "Taxa-Path-Spear_CORRELATION_COEFFICIENT", sep = "\t", row.names = FALSE, col.names =TRUE, append = FALSE)
  write.table(Corr_P_Value_Results, file = "Taxa-Path-Spear_CORRELATION_P_VALUE", sep = "\t", row.names = FALSE, col.names =TRUE, append = FALSE)
}