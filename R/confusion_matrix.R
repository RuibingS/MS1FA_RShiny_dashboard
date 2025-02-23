library(ggplot2)
library(dplyr)
library(RColorBrewer)



conf_matrix_plot_fun <- function(output_dir, input_list,title_txt) {
 
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  
  lapply(names(input_list), function(name) {
    
    x <- input_list[[name]] 

    if (!is.data.frame(x$Confusion_Matrix)) {
      stop(paste("Confusion_Matrix for", name, "is not a valid data frame."))
    }

    conf_matrix <- data.frame(
      Prediction = factor(rep(c("Positive", "Negative"), each = 2), 
                          levels = rev(c("Negative", "Positive"))),
      Actual = factor(rep(c("Positive", "Negative"), 2), 
                      levels = rev(c("Positive", "Negative"))),
      Count = as.numeric(as.vector(unlist(x$Confusion_Matrix)))  # Ensure numeric counts
    )

    palette_colors <- brewer.pal(9, "Blues")[4:7]

    p <- ggplot(conf_matrix, aes(x = Prediction, y = Actual, fill = Count)) +
      geom_tile(color = "white", linewidth = 1) + 
      geom_text(aes(label = Count), size = 6, fontface = "bold", color = "black") +
      scale_fill_gradientn(colors = palette_colors) +
      labs(title = paste0(title_txt, " ",name),
           x = "Predicted Label",
           y = "Actual Label") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "none",
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white", color = "white"), 
        plot.background = element_rect(fill = "white", color = "white")
      )
    output_filename <- file.path(output_dir, paste0(name, ".png"))
 
    
    #ggsave(output_filename, plot = p, width = 8, height = 8, units = "in", dpi = 300)
    
    ggsave(output_filename, plot = p, width = 5, height = 5, dpi = 300)
    
    return(p)  
  })
}


#  test 
# conf_matrix_plot_fun(output_dir=here::here("output_confusion_matrix") , conf_matrix_MS1FA_mzmine, title_txt=paste("Confusion Matrix - MS1FA"))
# 
# conf_matrix_plot_fun(output_dir=here::here("output_confusion_matrix"), conf_matrix_MS1FA_XCMS, title_txt=paste("Confusion Matrix - MS1FA"))
#  
# conf_matrix_plot_fun(output_dir=here::here("output_confusion_matrix"),conf_matrix_MZmine_results, title_txt=paste("Confusion Matrix - MZmine"))
#  
# conf_matrix_plot_fun(output_dir=here::here("output_confusion_matrix"),CAMERA_results, title_txt=paste("Confusion Matrix - CAMERA"))
# 

