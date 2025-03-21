height_interp <- function(height_est,h,parameter){
  interpoints = height_est$df
  interpoints = split(interpoints,interpoints$group)
  
  interpoints = simplify2array(interpoints)
  
  cn = colnames(interpoints)
  rn = rownames(interpoints)
  
  interpoints = array(unlist(interpoints),dim = c(length(interpoints[[1]]),length(rn),length(cn)))
  
  # Original data
  
  p <-  interpoints[,1,1]
  data <- interpoints[,2,]
  
  
  ############################## INTERPOLATION
  # New h values
  h_new <- seq(10, 300, 1) # Replace 100 with the number of new h values you want
  # Converting the matrix to vectors
  h_expanded <- rep(h, each = length(p))
  p_expanded <- rep(p, times = length(h))
  z_values <- as.vector(data)  # Transpose to match the (h, p) structure
  
  # Interpolating along the h dimension
  interpolation_result <- akima::interp(x = h_expanded, y = p_expanded, z = z_values, xo = h_new, yo = p)
  
  interpolated_data <- matrix(interpolation_result$z, nrow = length(h_new), ncol = length(p))
  
  ##############################
  
  
  
  # Plotting the original and interpolated data for visualization (optional)
  par(mfrow = c(1, 2))
  
  fields::image.plot(h, p, t(data), main = parameter, xlab = "height", ylab = "sample fraction", log = "y")
  
  # Plot the interpolated data with a color scale and log y-axis
fields::image.plot(h_new, p, interpolated_data, main = paste("Interpolated",parameter), xlab = "height", ylab = "sample farction", log = "y")
  return(interpolated_data)
  }

height_heatmap <- function(height_est,h,parameter){
  interpoints = height_est$df
  interpoints = split(interpoints,interpoints$group)
  
  interpoints = simplify2array(interpoints)
  
  cn = colnames(interpoints)
  rn = rownames(interpoints)

  interpoints = array(unlist(interpoints),dim = c(length(interpoints[[1]]),length(rn),length(cn)))
  

  
  
  zlim = range(interpoints[,c(2,3,4),])
  
  data = interpoints[,2,]
  p <-  interpoints[,1,1]
  rownames(data) <- rep("", length(p))
  
  rownames(data)[seq(1, length(p), by = 5)] <- round(p,digits = 3)[seq(1, length(p), by = 5)]
  
  colnames(data) <- h
  
  
  # Define the heatmap
  ht <- Heatmap(data,
                name = parameter,
                column_title = "height",
                row_title = "sample fraction",
                row_names_side = "left",
                column_names_side = "top",
                col = colorRamp2(c(zlim[1], mean(zlim), zlim[2]), c("blue", "yellow", "red")),
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_column_names = TRUE,
                show_row_names = TRUE,
                heatmap_legend_param = list(title = parameter))
  
  
  bp = cbind(interpoints[,3,],interpoints[,4,])
  
  # Create bar annotation
  ha <- rowAnnotation(
    confidence = anno_boxplot(bp, border = FALSE, gp = gpar(fill = colorRamp2(c(zlim[1], mean(zlim), zlim[2]), c("blue", "yellow", "red"))(bp)),height = unit(10,"cm")),
    annotation_name_side = "top"
  )
  
  
  # Draw the heatmap with annotation
  draw(ht+ha, heatmap_legend_side = "right", annotation_legend_side = "right")
}