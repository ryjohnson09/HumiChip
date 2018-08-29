####################################################
# Name: ordination_plot.R
# Date Created: 27 August 2018
# Author: Ryan Johnson
# Purpose: This scriptcontains a function 
#   that can be used to plot the results from the
#   ordination_function.R function (the ord_object)
####################################################

ord_plot <- function(ord_object, # Results from running `ordination_function.R`
                     ord_type, # PCoA, PCA, or DCA
                     color_points_by = visit_number,
                     connect_points_by_studyID = FALSE,
                     plot_title = "",
                     y_lab = "", 
                     x_lab = ""){
  # Libraries
  library(ggplot2)
  library(ggExtra)
  
  # Determine ord type for axis labels
  if (ord_type == "PCoA"){
    x_axis = "Axis.1"
    y_axis = "Axis.2"
  } else if (ord_type == "PCA"){
    x_axis = "PC1"
    y_axis = "PC2"
  } else if (ord_type == "DCA"){
    x_axis = "DCA1"
    y_axis = "DCA2"
  } else {
    stop("Invalid x/y axes")
  }
  
  
  # Plotting
  ord_plot <- ggplot(ord_object$humichip_data, 
                     aes_string(x = x_axis, y = y_axis)) +
    xlab(x_lab) +
    ylab(y_lab) +
    geom_point(aes_string(color = color_points_by), alpha = 0.8, size = 6) +
    geom_point(pch = 21, size = 6, alpha = 0.6) +
    ggtitle(plot_title) + 
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 14), 
      axis.text.x = element_text(size = 12, hjust = 1),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 14), 
      plot.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 12),
      legend.title = element_blank())
  
  # Connet points?
  if(connect_points_by_studyID == TRUE){
    ord_plot <- ord_plot + geom_line(aes(group = STUDY_ID), linetype = 1)
  }
  
  ggMarginal(ord_plot, groupFill = TRUE, groupColour = TRUE)
}