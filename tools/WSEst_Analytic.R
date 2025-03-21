########################################################################################

# WIND SPEED ESTIMATION WITH ANALYTICAL UNCERTAINTIES FOR ALL METHODS AND DISTRIBUTIONS
# METHODS CHOSEN FOR THE ANALYSIS

########################################################################################

# ==================================================
# GEV DISTRIBUTION
#
{
  
# Function to estimate wind speed using Block Maxima (BM) approach
# Compares model simulation data with measured data for tail index estimation
# Parameters:
#   data_model: The model simulation data, typically a data frame or matrix with wind speeds
#   data_measure: The measurement data, typically a data frame or matrix with wind speeds
#   col_model: The column name (or index) in data_model that contains the wind speed data (default: 'F010')
#   col_measure: The column name (or index) in data_measure that contains the wind speed data (default: 'F010')
#   tag: Optional parameter for tagging or labeling the output (default: NULL)
#   length.out: Number of points to compute for return periods (default: 30)
#   cov_type: The type of covariance matrix used for fitting ("num","analytic")
#   fixtailto0: Boolean flag to fix the tail index to zero. If TRUE, the tail index is set to zero (default: FALSE)

WSEst_model_to_measure_BM_4 <- function(data_model, data_measure, col_model='F010', col_measure='F010', tag=NULL, length.out=30, cov_type="num", fixtailto0=F) {
  
  # Conditionally set tail to zero or fit the model for tail index
  if (fixtailto0){
    tail_model <- list(tail = 0, tailStd = 0)
  } else {
    # Select block maxima from model data
    bm_model <- BM_select(data_model, height = col_model)
    
    # Fit the BM model to the data to get the tail index
    tail_model <- BM_fit_cov(bm_model, plot = F, cov_type=cov_type)
    
    # Extract the tail and its standard deviation
    tail_model <- list(tail = tail_model$trueTail['tail'], tailStd = tail_model$cov['tail', 'tail'])
  }
  
  # Select block maxima from the measurement data
  bm <- BM_select(data_measure, height = col_measure)
  
  # Function to constrain the tail index parameter during optimization
  fpar <- function(p, xpar) {
    loc <- matrix(p[1], xpar[1], xpar[2])
    scale <- matrix(p[2], xpar[1], xpar[2])
    shape <- matrix(tail_model$tail, xpar[1], xpar[2])
    list(loc = loc, scale = scale, shape = shape)
  }
  
  # Initial guesses for location and scale parameters
  kwargs <- list(start = c(mean(bm$max), sd(bm$max)), fpar = fpar, xpar = c(length(bm$max), 1), tailStd = tail_model$tailStd)
  
  # Quantile for 95% confidence intervals
  z <- qnorm(1 - 0.05 / 2)
  
  # Fit the measurement data with BM model (without using model's tail)
  fit_measure <- BM_fit_cov(bm, plot = F, cov_type=cov_type)
  
  # Fit the measurement data with BM model, including the model's tail index
  fit_measure_model <- BM_fit_cov(bm, plot = F, kwargs = kwargs, cov_type=cov_type)
  
  # Define the range of return periods to be evaluated
  X <- 10^seq(log10(5), log10(10000), length.out = length.out)
  X <- sort(unique(c(c(50, 100, 1000, 10000), X)))
  
  # Print model tail estimation for logging
  print(paste('Model estimation : tail =', tail_model$tail, 'sd', tail_model$tailStd))
  
  # Compute estimates based on original measurement data
  parameters = list(mean = fit_measure$trueTail, cov = fit_measure$cov)
  measure_est = qgev_distrib_2(exp(-1 / X), parameters)
  
  # Print original estimate parameters
  print('Original estimation : ')
  print(parameters)
  
  # Compute estimates based on the model's tail index and measurement data
  parameters = list(mean = fit_measure_model$trueTail, cov = fit_measure_model$cov)
  model_est = qgev_distrib_2(exp(-1 / X), parameters)
  
  # Print combined estimate (measurements + model tail)
  print('Measurements + model tail estimation : ')
  print(parameters)
  
  # Create a data frame for storing results (return periods and estimates)
  df <- data.frame(
    return = X,
    original = measure_est$mean,
    lb_o = measure_est$lb,
    ub_o = measure_est$ub,
    sd_o = measure_est$sd,
    model_est = model_est$mean,
    lb = model_est$lb,
    ub = model_est$ub,
    sd = model_est$sd,
    proba = 1 / X
  )
  
  # Create a data frame for observed points from the measurement data
  points <- data.frame(speed = sort(bm$max), return = 1 / (1 - seq(1 / length(bm$max), 1 - 1 / length(bm$max), length.out = length(bm$max))))
  
  # Plotting block maxima estimation using ggplot2
  {
    custom_color <- "blue"
    
    # Create the base plot for estimated values and confidence intervals
    gg <- ggplot(df, aes(x = return)) +
      geom_line(aes(y = model_est, color = 'shape estimated with model')) +
      geom_ribbon(aes(ymin = lb, ymax = ub, fill = 'Confidence Interval, 95%'), alpha = 0.2) +
      geom_line(aes(y = original, color = 'all parameters estimated\nfrom measurements')) +
      geom_line(aes(y = lb_o), lty = 'dashed', color = "red") +
      geom_line(aes(y = ub_o), lty = 'dashed', color = "red") +
      geom_point(data = points, aes(y = speed, shape = 'Measurements observed')) +
      scale_x_log10(labels = scales::comma) +
      ylim(range(c(df$original, df$model_est, df$lb, df$ub, points$speed))) +
      labs(x = 'Return Period (years)', y = 'Return Wind Speed (m/s)', title = paste0('Estimated wind speed for return period until 10000 years,\n1y BM (GEV), on measurements')) +
      scale_color_manual(values = c('shape estimated with model' = custom_color, 'all parameters estimated\nfrom measurements' = 'red', 'TRUE' = 'black', 'FALSE' = 'lightgrey')) +
      scale_fill_manual(values = c('Confidence Interval, 95%' = custom_color)) +
      theme(legend.position = "right")
    
    # Extract the last row of the data frame for annotation purposes
    last_values <- tail(df, 1)
    
    # Annotate the plot with the final values for original and model estimates
    gg <- gg +
      geom_text(data = last_values, aes(label = sprintf("%.2f", original), x = max(df$return), y = original), color = 'red', size = 3, hjust = 0) +
      geom_text(data = last_values, aes(label = sprintf("%.2f", model_est), x = max(df$return), y = model_est), color = custom_color, size = 3, hjust = 0) +
      geom_text(data = last_values, aes(label = sprintf("%.2f", lb), x = max(df$return), y = lb), color = custom_color, size = 3, hjust = 0) +
      geom_text(data = last_values, aes(label = sprintf("%.2f", ub), x = max(df$return), y = ub), color = custom_color, size = 3, hjust = 0)
    
    # Find the row closest to a return period of 50 years
    closest_to_50 <- df %>% slice(which.min(abs(return - 50)))
    
    # Highlight the values closest to a return period of 50 years
    gg <- gg +
      geom_point(data = closest_to_50, aes(x = return, y = original), color = "red", size = 1) +
      geom_point(data = closest_to_50, aes(x = return, y = model_est), color = "blue", size = 1) +
      geom_point(data = closest_to_50, aes(x = return, y = lb), color = "blue", size = 1) +
      geom_point(data = closest_to_50, aes(x = return, y = ub), color = "blue", size = 1) +
      geom_text(data = closest_to_50, aes(label = sprintf("%.2f", original), x = return, y = original + 1), color = "red", size = 3, hjust = 0) +
      geom_text(data = closest_to_50, aes(label = sprintf("%.2f", model_est), x = return, y = model_est + 1), color = "blue", size = 3, hjust = 0) +
      geom_text(data = closest_to_50, aes(label = sprintf("%.2f", lb), x = return, y = lb - 1), color = "blue", size = 3, hjust = 0) +
      geom_text(data = closest_to_50, aes(label = sprintf("%.2f", ub), x = return, y = ub + 1), color = "blue", size = 3, hjust = 0)
  }
  
  # Return a list containing the plot, data frame, observed points, and parameter distributions
  parameter_distributions = list(
    model = tail_model,
    measure = fit_measure,
    measure_model = fit_measure_model
  )
  return(list(gg = gg, df = df, timestep = 1, points = points, parameter_distributions = parameter_distributions))
}

}


# ==================================================
# GP DISTRIBUTION
#
{
  
# Function to estimate wind speed extremes using the Generalized Pareto (GP) distribution
# Parameters:
#   data: The dataset, typically containing wind speed and datetime columns
#   timestep: The timestep (in hours) used for the analysis
#   threshold: The percentage threshold to define the peaks over threshold (PoT)
#   col: The column name (or index) in 'data' that contains the wind speed data (default: 'F010')
#   length.out: The number of return period points to compute (default: 30)
#   fixtailto0: Boolean flag to fix the tail index to zero. If TRUE, the tail index is set to zero (GUMBEL) (default: FALSE)
#   winter: Boolean flag to indicate whether to consider only winter months (October-March) for the analysis (default: FALSE)

WSEst_GP_4 <- function(data, timestep, threshold, col = "F010", length.out = 30, fixtailto0 = F, winter = F) {
  
  # Linearly interpolate missing wind speed data in the specified column
  ws <- na.approx(data[[col]])
  
  # If 'winter' flag is TRUE, only use data from the winter months (October-March)
  if (winter) {
    ws <- na.approx(data[[col]][month(data$DateTime) %in% c(10:12, 1:3)])
  }
  
  # Select peaks over the threshold (PoT) with a threshold of 0.3 and minimum time between peaks of 12 hours (adjusted by timestep)
  pot <- PoTselect_2(ws, 0.3, 12 / timestep)
  
  # Print the number of peaks identified
  print(length(pot$pot))
  
  # Calculate the number of peaks to be used in the estimation (l0) based on the given threshold
  l0 <- round(threshold * length(pot$pot)) - 1
  
  # If 'fixtailto0' is TRUE, fix the tail index to 0, otherwise leave it NULL
  if (fixtailto0) {
    fixedpar <- list(gamma0 = 0, gamma0Std = 0)
  } else {
    fixedpar <- NULL
  }
  
  # Calculate the average time between peaks in years
  timestep_selected_peaks <- mean(diff(pot$ind) * timestep / (24 * 365.2425))
  
  # Create a data frame with sorted peak wind speeds and their corresponding probabilities
  points <- data.frame(speed = sort(pot$pot))
  points$p <- seq(0, 1 - 1 / length(points$speed), length.out = length(points$speed))  # Calculate empirical probabilities
  points$used <- 1 - points$p <= l0 / length(pot$pot)  # Mark peaks that will be used in the model
  points$return <- timestep_selected_peaks / (1 - points$p)  # Calculate return periods for the observed peaks
  
  # Create return periods for which the model will predict wind speeds
  X <- 10^seq(log10(timestep_selected_peaks * length(pot$pot) / l0), log10(10000), length.out = length.out)
  X <- sort(unique(c(c(50, 100, 1000, 10000), X)))  # Ensure standard return periods are included
  
  # Fit the Generalized Pareto distribution to the peaks using maximum likelihood estimation
  pot_fit <- FitGP_MLE2(pot$pot, timestep_selected_peaks / X, N = 0, r11 = 1, fixedpar = fixedpar, l0 = l0, sigma = Inf)
  
  # Create a data frame of the estimated mean wind speeds and confidence intervals for each return period
  distrib <- data.frame(
    mean = t(pot_fit$quantile),
    lb = t(pot_fit$quantile - qnorm(0.975) * pot_fit$quantileStd),  # Lower bound of 95% confidence interval
    ub = t(pot_fit$quantile + qnorm(0.975) * pot_fit$quantileStd),  # Upper bound of 95% confidence interval
    sd = t(pot_fit$quantileStd)
  )
  
  # Prepare the data frame for plotting return periods and wind speeds
  df <- data.frame(
    return = X,
    original = distrib$mean,
    lb = distrib$lb,
    ub = distrib$ub,
    sd = distrib$sd
  )
  
  # Plot the estimated wind speeds and confidence intervals
  gg <- ggplot(df, aes(x = return)) +
    geom_line(aes(y = original), col = 'red') +  # Line for estimated wind speeds
    geom_ribbon(aes(ymin = lb, ymax = ub), col = 'red', alpha = 0.2) +  # Shaded confidence interval
    geom_point(data = points, aes(x = return, y = speed, color = used)) +  # Observed peaks
    geom_vline(xintercept = timestep_selected_peaks * length(pot$pot) / l0, lty = 'dashed', alpha = 0.3) +  # Vertical dashed line for the threshold
    scale_x_log10(labels = scales::comma) +  # Logarithmic scale for return periods
    labs(x = 'Return Period (years)', y = 'Wind Speed (m/s)',
         title = paste0('Estimated wind speed for return period until 10000 years,\nPoT: ', threshold * 100, '% of the peaks, GP')) +
    theme(legend.position = "right") +
    scale_color_manual(name = "", values = c('TRUE' = 'black', 'FALSE' = 'lightgrey'))  # Color the used/unused points
  
  # Get the last values of the return period estimates for annotation
  last_values <- tail(df, 1)
  print(last_values)
  
  # Annotate the last values (largest return period) on the plot
  gg <- gg +
    geom_text(data = last_values, aes(label = sprintf("%.2f", original), x = max(df$return), y = original), 
              color = 'red', size = 3, hjust = 0) +
    geom_text(data = last_values, aes(label = sprintf("%.2f", lb), x = max(df$return), y = lb), 
              color = 'red', size = 3, hjust = 0) +
    geom_text(data = last_values, aes(label = sprintf("%.2f", ub), x = max(df$return), y = ub), 
              color = 'red', size = 3, hjust = 0)
  
  # Find the closest return period to 50 years and highlight it on the plot
  closest_to_50 <- df %>%
    slice(which.min(abs(return - 50)))
  
  # Annotate the values near the 50-year return period
  gg <- gg +  
    geom_point(data = closest_to_50, aes(x = return, y = original), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = lb), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = ub), color = "blue", size = 1) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", original), x = return, y = original + 1), 
              color = "blue", size = 3, hjust = 0) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", lb), x = return, y = lb - 1), 
              color = "blue", size = 3, hjust = 0) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", ub), x = return, y = ub + 1), 
              color = "blue", size = 3, hjust = 0)
  
  # Return the plot, data, and points used in the analysis
  return(list(gg = gg, df = df, points = points))
}




# Function to estimate wind speed using Peak over Threshold (POT) approach with Generalized Pareto (GP) distribution
# Compare tail index estimation on model simulations and measured wind speeds using Generalized Pareto distribution (GP)
# Parameters:
#   data_model: Dataset from the model simulations containing wind speed data and datetime columns
#   data_measure: Dataset of measured wind speed values with datetime columns
#   col_model: The column in data_model representing wind speed (default: 'PF010')
#   col_measure: The column in data_measure representing wind speed (default: 'F010')
#   timestep_model: The timestep of the model data (in hours)
#   timestep_measure: The timestep of the measured data (in hours)
#   th_model: Threshold for the model to define peaks over threshold (PoT)
#   th_measure: Threshold for the measured data to define PoT
#   length.out: The number of return period points to compute (default: 100)
#   peak_frac: Percentage of the storm threshold used in PoT selection (default: 0.3)
#   winter: Boolean flag to consider only winter months (default: FALSE)

WSEst_model_to_measure_GP_4 <- function(data_model, data_measure, col_model='PF010', col_measure='F010', timestep_model=1, timestep_measure=1, th_model=NULL, th_measure=NULL, length.out=100, peak_frac=0.3, winter=F) {
  
  # Linearly interpolate missing wind speed data in the model and measurement datasets
  ws_model <- na.approx(data_model[[col_model]])
  ws_measure <- na.approx(data_measure[[col_measure]])
  
  # If 'winter' flag is TRUE, select only winter months (October-March)
  if (winter) {
    ws_model <- na.approx(data_model[[col_model]][month(data_model$DateTime) %in% c(10:12, 1:3)])
    ws_measure <- na.approx(data_measure[[col_measure]][month(data_measure$DateTime) %in% c(10:12, 1:3)])
  }
  
  # PoT selection for the model
  cat("Estimation on model\n")
  flush.console()
  pot_model <- PoTselect_2(ws_model, peak_frac, 12/timestep_model)
  
  # Calculate the number of peaks to use based on the threshold for the model
  l0_model <- round(th_model * length(pot_model$pot)) - 1
  
  # Fit the GP distribution to the model data using maximum likelihood estimation (MLE)
  fit_model <- FitGP_MLE2(pot_model$pot, 1, N=nrow(data_model), r11=1, l0=l0_model)
  
  # Fix the tail index for measurements to match the tail index estimated from the model
  fixedpar <- list(gamma0=fit_model$tailindex, gamma0Std=fit_model$tailindexStd)
  
  # PoT selection for the measured data
  pot_measure <- PoTselect_2(ws_measure, peak_frac, 12/timestep_measure)
  print(sum(is.na(pot_measure$pot)))
  
  # Calculate the number of peaks to use based on the threshold for measurements
  l0_measure <- round(th_measure * length(pot_measure$pot)) - 1
  
  # Calculate the average time between peaks in years for the measured data
  timestep_POT_measure <- (length(data_measure[[col_measure]]) / length(pot_measure$pot)) * timestep_measure / (24 * 365.2425)
  
  # Define return periods to estimate wind speeds
  X <- 10^seq(log10(timestep_POT_measure * length(pot_measure$pot) / l0_measure), log10(10000), length.out=length.out)
  X <- sort(unique(c(c(50, 100, 1000, 10000), X)))  # Include standard return periods
  
  # Fit GP distribution to the measured data
  fit_measure <- FitGP_MLE2(pot_measure$pot, p=timestep_POT_measure/X, r11=1, l0=l0_measure)
  
  # Fit GP distribution to measured data, but fix the tail index to the model's estimate
  fit_measure_model <- FitGP_MLE2(pot_measure$pot, 1, N=nrow(data_measure), fixedpar=fixedpar, r11=1, l0=l0_measure)
  
  # Print the model's GP parameters
  print(paste('Model estimation: tail =', fit_model$tailindex, 'sd', fit_model$tailindexStd, 
              'scale =', fit_model$scale, 'loc =', fit_model$location))
  
  # Print the GP parameters estimated from measurements without the fixed tail index
  loc <- list(mean=fit_measure$location, sd=fit_measure$locationStd)
  scale <- list(mean=fit_measure$scale, sd=fit_measure$logdispStd)
  tail <- list(mean=fit_measure$tailindex, sd=fit_measure$tailindexStd)
  print(paste('Original estimation (without fixed tail): tail =', tail$mean, 'sd', tail$sd, 'scale =', scale$mean, 'loc =', loc$mean))
  
  # Create a data frame with estimated quantiles and confidence intervals from measured data
  measure_est <- data.frame(
    mean=t(fit_measure$quantile),
    lb=t(fit_measure$quantile - qnorm(0.975) * fit_measure$quantileStd),
    ub=t(fit_measure$quantile + qnorm(0.975) * fit_measure$quantileStd),
    sd=t(fit_measure$quantileStd)
  )
  
  # Print the GP parameters when using the fixed tail index from the model
  loc <- list(mean=fit_measure_model$location, sd=fit_measure_model$locationStd)
  scale <- list(mean=fit_measure_model$scale, sd=fit_measure_model$logdispStd)
  tail <- list(mean=fit_measure_model$tailindex, sd=fit_measure_model$tailindexStd)
  print(paste('Measurements + Model estimation: tail =', tail$mean, 'sd', tail$sd, 'scale =', scale$mean, 'loc =', loc$mean))
  
  # Estimate return values using the model-estimated tail index
  model_est <- qgp_distrib(p=1 - timestep_POT_measure / X, loc=loc, scale=scale, tail=tail, t=length(pot_measure$pot) / l0_measure)
  
  # Create a data frame for plotting with return periods, estimated wind speeds, and confidence intervals
  df <- data.frame(
    return=X,
    original=measure_est$mean,
    lb_o=measure_est$lb,
    ub_o=measure_est$ub,
    sd_o=measure_est$sd,
    model_est=model_est$mean,
    lb=model_est$lb,
    ub=model_est$ub,
    sd=model_est$sd,
    proba=timestep_POT_measure / X
  )
  
  # Create a data frame for observed peak values
  points <- data.frame(speed=sort(pot_measure$pot))
  points$p <- seq(0, 1 - 1/length(points$speed), length.out=length(points$speed))  # Observed probabilities
  points$used <- 1 - points$p <= l0_measure / length(pot_measure$pot)  # Mark peaks used for estimation
  points$return <- timestep_POT_measure / (1 - points$p)  # Return periods for observed peaks
  
  # Create the ggplot object for visualizing return periods and wind speeds
  custom_color <- "blue"
  gg <- ggplot(df, aes(x=return)) +
    geom_line(aes(y=model_est, color='shape estimated with model')) +
    geom_ribbon(aes(ymin=lb, ymax=ub, fill='Confidence Interval, 95%'), alpha=0.2) +
    geom_line(aes(y=original, color='all parameters estimated\nfrom measurements')) +
    geom_line(aes(y=lb_o), lty='dashed', color="red") +
    geom_line(aes(y=ub_o), lty='dashed', color="red") +
    geom_point(data=points, aes(y=speed, color=used, shape='Measurements observed')) +
    scale_x_log10(labels=scales::comma) +
    ylim(range(c(df$original, df$model_est, df$lb, df$ub, points$speed))) +
    labs(x='Return Period (years)', y='Return Wind Speed (m/s)', 
         title=paste0('Estimated wind speed for return period until 10000 years,\nPoT (GP), p=', peak_frac, ', th=', th_measure, ', on measurements')) +
    scale_color_manual(values=c('shape estimated with model'=custom_color, 'all parameters estimated\nfrom measurements'='red', 'TRUE'='black', 'FALSE'='lightgrey')) +
    scale_fill_manual(values=c('Confidence Interval, 95%'=custom_color)) +
    theme(legend.position="right")
  
  # Annotate the last values on the plot
  last_values <- tail(df, 1)
  print(last_values)
  gg <- gg +
    geom_text(data=last_values, aes(label=sprintf("%.2f", original), x=max(df$return), y=original), color='red', size=3, hjust=0) +
    geom_text(data=last_values, aes(label=sprintf("%.2f", model_est), x=max(df$return), y=model_est), color=custom_color, size=3, hjust=0) +
    geom_text(data=last_values, aes(label=sprintf("%.2f", lb), x=max(df$return), y=lb), color=custom_color, size=3, hjust=0) +
    geom_text(data=last_values, aes(label=sprintf("%.2f", ub), x=max(df$return), y=ub), color=custom_color, size=3, hjust=0)
  
  # Highlight the values closest to the 50-year return period
  closest_to_50 <- df %>% slice(which.min(abs(return - 50)))
  
  gg <- gg +  
    geom_point(data=closest_to_50, aes(x=return, y=original), color="red", size=1) +
    geom_point(data=closest_to_50, aes(x=return, y=model_est), color="blue", size=1) +
    geom_point(data=closest_to_50, aes(x=return, y=lb), color="blue", size=1) +
    geom_point(data=closest_to_50, aes(x=return, y=ub), color="blue", size=1) +
    geom_text(data=closest_to_50, aes(label=sprintf("%.2f", original), x=return, y=original + 1), color="red", size=3, hjust=0) +
    geom_text(data=closest_to_50, aes(label=sprintf("%.2f", model_est), x=return, y=model_est + 1), color="blue", size=3, hjust=0) +
    geom_text(data=closest_to_50, aes(label=sprintf("%.2f", lb), x=return, y=lb - 1), color="blue", size=3, hjust=0) +
    geom_text(data=closest_to_50, aes(label=sprintf("%.2f", ub), x=return, y=ub + 1), color="blue", size=3, hjust=0)
  
  # Return the plot, data frame, timestep, observed points, and parameter distributions
  parameter_distributions <- list(
    model=list(tail=list(mean=fit_model$tailindex, sd=fit_model$tailindexStd),
               scale=list(mean=fit_model$scale, sd=fit_model$logdispStd),
               loc=list(mean=fit_model$location, sd=fit_model$locationStd)),
    measure=fit_measure,
    measure_model=list(tail=tail, scale=scale, loc=loc)
  )
  
  return(list(gg=gg, df=df, timestep=timestep_POT_measure, points=points, parameter_distributions=parameter_distributions))
}

}


# ==================================================
# GW DISTRIBUTION
#
{

# wind speed estimation for POT method
#   fixtailto1: Boolean flag to fix the tail index to one. If TRUE, the tail index is set to 1 (GUMBEL) (default: FALSE)
  
WSEst_GW_2 <- function(data,timestep,threshold,col="F010",length.out = 30,fixtailto1=F,winter=F){
  ws <- na.approx(data[[col]])
  if (winter){
    ws <- na.approx(data[[col]][month(data$DateTime) %in% c(10:12,1:3)])
  }
  pot <- PoTselect_2(ws,0.3,12/timestep)
  l0 <- round(threshold*length(pot$pot))-1
  
  if (fixtailto1){fixedpar <- list(theta0 = 1,theta0Std = 0)}
  else{fixedpar <- NULL}
  
  pot_fit <- FitGW_iHilli(pot$pot,p=1,l0 = l0,fixedpar = fixedpar)
  
  

  timestep_selected_peaks = mean(diff(pot$ind)*timestep/(24*365.2425)) # mean of the timediff between each peak
  
  
  points <- data.frame(speed = sort(pot$pot))
  points$p <- seq(0,1-1/length(points$speed),length.out=length(points$speed))     # observed Values
  points$used <- 1-points$p <= l0/length(pot$pot)
  points$return <- timestep_selected_peaks / (1- points$p)
  
  
  # print(ggplot(points,aes(x=speed,y=p,color=used))+geom_line()+geom_point()+scale_color_manual(name="",values = c('TRUE'='black','FALSE'='lightgrey')))
  
  
  tail <- list(mean=pot_fit$tailindex,sd=pot_fit$tailindexStd)
  print(tail)
  scale <- list(mean=pot_fit$scale,sd=pot_fit$logdispStd)
  loc <- list(mean=pot_fit$location,sd=pot_fit$locationStd)
  
  X <- 10^seq(log10(timestep_selected_peaks*length(pot$pot)/l0),log10(10000),length.out = length.out)
  X <- sort(unique(c(c(50,100,1000,10000),X)))
  
  
  distrib = qgw_distrib(1-timestep_selected_peaks/X,tail = tail,scale = scale,loc = loc,y=pot_fit$y)
  
  
  df <- data.frame(
    return = X,
    original = distrib$mean,
    lb = distrib$lb,
    ub = distrib$ub,
    sd = distrib$sd
  )
  
  
  gg = ggplot(df, aes(x = return)) +
    geom_line(aes(x=return,y=original),col='red') +
    geom_ribbon(aes(ymin=lb,ymax=ub),col='red',alpha=0.2)+
    geom_point(data=points,aes(x = return,y = speed,color=used)) +
    geom_vline(xintercept = timestep_selected_peaks*length(pot$pot)/l0,lty='dashed',alpha=0.3)+
    scale_x_log10(labels= scales::comma) +
    labs(x = 'Return Period (years)',
         y = 'Wind Speed (m/s)',
         title = paste0('Estimated wind speed for return period until 10000 years,\nPoT: ',threshold*100,'% of the peaks, GW')) +
    theme(legend.position = "right")+
    scale_color_manual(name="",values = c('TRUE'='black','FALSE'='lightgrey'))
  
  last_values <- tail(df, 1)  # Get the last row of df_measure
  print(last_values)
  # Annotate the last values
  gg <- gg +
    geom_text(data = last_values, aes(label = sprintf("%.2f", original), 
                                      x = max(df$return), y = original ), 
              color = 'red', size = 3, hjust = 0) +
    geom_text(data = last_values, aes(label = sprintf("%.2f", lb), 
                                      x = max(df$return), y = lb), 
              color = 'red', size = 3, hjust = 0) +
    geom_text(data = last_values, aes(label = sprintf("%.2f", ub), 
                                      x = max(df$return), y = ub), 
              color = 'red', size = 3, hjust = 0)
  
  closest_to_50 <- df %>% 
    slice(which.min(abs(return - 50)))
  
  gg <- gg +  
    geom_point(data = closest_to_50, aes(x = return, y = original), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = lb), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = ub), color = "blue", size = 1) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", original), 
                                        x = return, y = original + 1), 
              color = "blue", size = 3, hjust = 0) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", lb), 
                                        x = return, y = lb - 1), 
              color = "blue", size = 3, hjust = 0) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", ub), 
                                        x = return, y = ub + 1), 
              color = "blue", size = 3, hjust = 0)
  return(list(gg=gg,df=df,points=points))
}




# Function to estimate wind speed using Peak over Threshold (POT) approach with Generalized Weibull (GW) distribution
# Compare tail index estimation on model simulations and measured wind speeds using Generalized Weibull distribution (GW)
# Parameters:
#   data_model: Dataset from the model simulations containing wind speed data and datetime columns
#   data_measure: Dataset of measured wind speed values with datetime columns
#   col_model: The column in data_model representing wind speed (default: 'PF010')
#   col_measure: The column in data_measure representing wind speed (default: 'F010')
#   timestep_model: The timestep of the model data (in hours)
#   timestep_measure: The timestep of the measured data (in hours)
#   th_model: Threshold for the model to define peaks over threshold (PoT)
#   th_measure: Threshold for the measured data to define PoT
#   length.out: The number of return period points to compute (default: 100)
#   peak_frac: Percentage of the storm threshold used in PoT selection (default: 0.3)
#   winter: Boolean flag to consider only winter months (default: FALSE)

WSEst_model_to_measure_GW_2 <- function(data_model, data_measure, col_model='PF010', col_measure='F010', timestep_model=1, timestep_measure=1, th_model=NULL, th_measure=NULL,length.out = 100, tag = NULL, peak_frac=0.3, winter=F) {
  
  ws_model <- na.approx(data_model[[col_model]])
  ws_measure <- na.approx(data_measure[[col_measure]])
  if (winter){
    ws_model <- na.approx(data_model[[col_model]][month(data_model$DateTime) %in% c(10:12,1:3)])
    ws_measure <- na.approx(data_measure[[col_measure]][month(data_measure$DateTime) %in% c(10:12,1:3)])
  }
  
  # Bootstrap on model to get uncertainties of tail index estimation
  cat("Estimation on model\n")
  flush.console()
  pot_model <- PoTselect_2(ws_model,peak_frac,12/timestep_model)
  l0_model <- round(th_model*length(pot_model$pot))-1
  fit_model <- FitGW_iHilli(pot_model$pot,p= 1, l0= l0_model)
  
  
  
  
  # Function to constrain the tail index parameter in the estimation
  fixedpar = list(theta0 = fit_model$tailindex,theta0Std = fit_model$tailindexStd)
  
  # Compute scale and loc for original and tail_index model estimated
  
  
  pot_measure <- PoTselect_2(ws_measure,peak_frac,12/timestep_measure)
  l0_measure <- round(th_measure*length(pot_measure$pot))-1
  fit_measure <- FitGW_iHilli(pot_measure$pot,p= 1, l0= l0_measure)
  fit_measure_model <- FitGW_iHilli(pot_measure$pot,p= 1,fixedpar = fixedpar, l0= l0_measure)
  
  
  
  # mean timestep for POT select on measurements
  timestep_POT_measure <- (length(data_measure[[col_measure]])/length(pot_measure$pot))*timestep_measure/(24*365.2425)
  
  # Return periods
  X <- 10^seq(log10(timestep_POT_measure*length(pot_measure$pot)/l0_measure), log10(10000), length.out =length.out )
  X <- sort(unique(c(c(50,100,1000,10000),X)))
  
  print(paste('Model estimation : tail =',fit_model$tailindex,'sd',fit_model$tailindexStd,'scale =',fit_model$scale,' loc =',fit_model$location))
  # measurement original estimation
  loc = list(mean = fit_measure$location, sd = fit_measure$locationStd)
  scale = list(mean = fit_measure$scale, sd = fit_measure$logdispStd)
  tail = list(mean = fit_measure$tailindex, sd = fit_measure$tailindexStd)
  
  
  print(paste('Original estimation : tail =',tail$mean,'sd',tail$sd,'scale =',scale$mean,' loc =',loc$mean))
  
  measure_est = qgw_distrib(1-timestep_POT_measure/X,tail = tail,scale = scale,loc = loc, y= fit_measure$y)
  
  fit_measure = list(tail=tail,scale=scale,loc=loc,y=fit_measure$y)
  
  
  # Estimate return values using the model-estimated tail index
  loc = list(mean = fit_measure_model$location, sd = fit_measure_model$locationStd)
  scale = list(mean = fit_measure_model$scale, sd = fit_measure_model$logdispStd)
  tail = list(mean = fit_measure_model$tailindex, sd = fit_measure_model$tailindexStd)
  print(paste('Measurements + Model estimation : tail =',tail$mean,'sd',tail$sd,'scale =',scale$mean,' loc =',loc$mean))
  
  model_est <- qgw_distrib(p = 1 - timestep_POT_measure / X, loc = loc, scale = scale, tail = tail, y= fit_measure_model$y)
  
  
  
  # Create a data frame for plotting
  df <- data.frame(
    return = X,
    original = measure_est$mean,
    lb_o = measure_est$lb,
    ub_o = measure_est$ub,
    sd_o = measure_est$sd,
    model_est = model_est$mean,
    lb = model_est$lb,
    ub = model_est$ub,
    sd = model_est$sd,
    proba = timestep_POT_measure/X
  )
  
  
  # Create a data frame for observed values
  
  points <- data.frame(speed = sort(pot_measure$pot))
  points$p <- seq(0,1-1/length(points$speed),length.out=length(points$speed))     # observed Values
  points$used <- 1-points$p <= l0_measure/length(pot_measure$pot)
  points$return <-  timestep_POT_measure / (1- points$p)
  
  
  
  
  # Create the ggplot
  custom_color <- "blue"
  
  gg <- ggplot(df, aes(x = return)) +
    geom_line(aes(y = model_est, color = 'shape estimated with model')) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = 'Confidence Interval, 95%'), alpha = 0.2) +
    geom_point(data = points, aes(y = speed, color=used,shape = 'Measurements observed')) +
    geom_line(aes(y = original, color = 'all parameters estimated\nfrom measurements')) +
    geom_line(aes(y = lb_o), lty = 'dashed', color = "red") +
    geom_line(aes(y = ub_o), lty = 'dashed', color = "red") +
    scale_x_log10(labels = scales::comma) +
    ylim(range(c(df$original, df$model_est, df$lb, df$ub, points$speed))) +
    labs(x = 'Return Period (years)', y = 'Return Wind Speed (m/s)', title = paste0('Estimated wind speed for return period until 10000 years,\nPoT (GW), p=',peak_frac,', th=',th_measure,', on measurements')) +
    scale_color_manual(values = c('shape estimated with model' = custom_color, 'all parameters estimated\nfrom measurements' = 'red','TRUE'='black','FALSE'='lightgrey')) +
    scale_fill_manual(values = c('Confidence Interval, 95%' = custom_color)) +
    theme(legend.position = "right")
  
  last_values <- tail(df, 1)  # Get the last row of df_measure
  # print(last_values)
  # Annotate the last values
  gg <- gg +
    geom_text(data = last_values, aes(label = sprintf("%.2f", original), 
                                      x = max(df$return), y = original ), 
              color = 'red', size = 3, hjust = 0) +
    geom_text(data = last_values, aes(label = sprintf("%.2f", model_est), 
                                      x = max(df$return), y = model_est), 
              color = custom_color, size = 3, hjust = 0) +
    geom_text(data = last_values, aes(label = sprintf("%.2f", lb), 
                                      x = max(df$return), y = lb), 
              color = custom_color, size = 3, hjust = 0) +
    geom_text(data = last_values, aes(label = sprintf("%.2f", ub), 
                                      x = max(df$return), y = ub), 
              color = custom_color, size = 3, hjust = 0)
  
  closest_to_50 <- df %>% 
    slice(which.min(abs(return - 50)))
  
  gg <- gg +  
    geom_point(data = closest_to_50, aes(x = return, y = original), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = model_est), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = lb), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = ub), color = "blue", size = 1) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", original), 
                                        x = return, y = original + 1), 
              color = "blue", size = 3, hjust = 0) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", model_est), 
                                        x = return, y = model_est + 1), 
              color = "blue", size = 3, hjust = 0) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", lb), 
                                        x = return, y = lb - 1), 
              color = "blue", size = 3, hjust = 0) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", ub), 
                                        x = return, y = ub + 1), 
              color = "blue", size = 3, hjust = 0)
  
  
  # Return a list containing the ggplot object, data frame, and observed points
  
  parameter_distributions = list(
    model = list(tail = list(mean = fit_model$tailindex,sd= fit_model$tailindexStd),
                 scale = list(mean = fit_model$scale,sd= fit_model$logdispStd),
                 loc = list(mean = fit_model$location,sd= fit_model$locationStd),
                 y=fit_model$y),
    measure = fit_measure,
    measure_model = list(tail=tail,scale=scale,loc=loc,y=fit_measure_model$y)
  )
  return(list(gg = gg, df = df,timestep = timestep_POT_measure,points = points, parameter_distributions = parameter_distributions))
}

}