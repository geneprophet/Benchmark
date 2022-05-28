default_labeller <- function(n) {
  function(str){
    str <- gsub("_", " ", str)
    yulab.utils::str_wrap(str, n)
  }
}

library(ggplot2)
dotplot = function(df,x="Dataset",y="Biomarker",size="AUC",color="p_value"){
  label_func <- default_labeller(40)
  p <- ggplot(df, aes_string(x = x, y = y, size = size)) +
    scale_y_discrete(labels = label_func)
  p <- p +  geom_point(aes_string(color = color)) 
  p + scale_color_continuous(low="red", high="blue",limits = c(0,0.05), breaks = c(0.00, 0.01, 0.02,0.03,0.04, 0.05),
                             guide=guide_colorbar(reverse=TRUE)) +
    ylab("Biomarker") + 
    # ggtitle("Benchmark") +
    DOSE::theme_dose(12) +
    theme(axis.text.x = element_text(angle = 90,hjust = 1)) +
    scale_size_continuous(range=c(3, 8)) + 
    guides(size  = guide_legend(order = 1), 
           color = guide_colorbar(order = 2))
  
}

