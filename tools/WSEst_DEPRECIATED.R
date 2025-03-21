
##################################################################################

# DEPRECIATED VERSIONS OF WIND SPEED EXTRAPOLATION

##################################################################################

# wind speed estimation for POT method with the first POT selection (deprecated, use version 1 of potselect which does not use the max of each excursion above the threshold)
# timestep in hours
# threshold : the percentage of data w want to use for the estimation
WSEst_GP <- function(data,timestep,threshold,col = "F010"){
  
  pot <- PoTselect(na.approx(data[[col]]),0.1,12/timestep)
  print(length(pot$pot))
  l0 <- round(threshold*length(pot$pot))-1
  pot_fit <- FitGP_MLE2(pot$pot,1, N= 0, r11= 1, fixedpar= NULL, l0= l0, sigma= Inf)
  
  # plot_extremes(pot,data)
  # plot(ecdf(pot$pot))
  # lines(seq(min(pot$pot),max(pot$pot),length.out=100),cgp((seq(min(pot$pot),max(pot$pot),length.out=100)-pot_fit$location)/pot_fit$scale,gamma = pot_fit$tailindex),col='red')
  
  
  timestep_selected_peaks = mean(diff(pot$ind)*timestep/(24*365.2425)) # mean of the timediff between each peak
  
  
  points <- data.frame(speed = sort(pot$pot))
  points$p <- seq(0,1-1/length(points$speed),length.out=length(points$speed))     # observed Values
  points$used <- 1-points$p <= l0/length(pot$pot)
  points$return <- timestep_selected_peaks / (1- points$p)
  
  
  # print(ggplot(points,aes(x=speed,y=p,color=used))+geom_line()+geom_point()+scale_color_manual(name="",values = c('TRUE'='black','FALSE'='lightgrey')))
  
  
  tail <- list(mean=pot_fit$tailindex,sd=pot_fit$tailindexStd)
  scale <- list(mean=pot_fit$scale,sd=pot_fit$logdispStd)
  loc <- list(mean=pot_fit$location,sd=pot_fit$locationStd)
  
  X <- 10^seq(log10(timestep_selected_peaks*length(pot$pot)/l0),log10(10000),length.out = 30)
  
  
  
  distrib = qgp_distrib(1-timestep_selected_peaks/X,tail = tail,scale = scale,loc = loc,t=length(pot$pot)/l0)
  
  
  df <- data.frame(
    return = X,
    original = distrib$mean,
    lb = distrib$lb,
    ub = distrib$ub
  )
  
  
  gg = ggplot(df, aes(x = return)) +
    geom_line(aes(x=return,y=original),col='red') +
    geom_ribbon(aes(ymin=lb,ymax=ub),col='red',alpha=0.2)+
    geom_point(data=points,aes(x = return,y = speed,color=used)) +
    geom_vline(xintercept = timestep_selected_peaks*length(pot$pot)/l0,lty='dashed',alpha=0.3)+
    scale_x_log10(labels= scales::comma) +
    labs(x = 'Return Period (years)',
         y = 'Wind Speed (m/s)',
         title = paste0('Estimated wind speed for return period until 10000 years,\nPoT: ',threshold*100,'% of the peaks, GP')) +
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





WSEst_model_to_measure_GP <- function(data_model, data_measure, col_model='PF010',col_measure='F010',timestep_model=1,timestep_measure=1, th_model=NULL, th_measure=NULL, tag=NULL) {
  
  # Bootstrap on model to get uncertainties of tail index estimation
  cat("Estimation on model\n")
  flush.console()
  pot_model <- PoTselect(data_model[[col_model]],0.1,12/timestep_model)
  l0_model <- round(th_model*length(pot_model$pot))-1
  fit_model <- FitGP_MLE2(pot_model$pot, 1, N= nrow(data_model), r11= 1, l0= l0_model)
  
  
  
  
  # Function to constrain the tail index parameter in the estimation
  fixedpar = list(gamma0 = fit_model$tailindex,gamma0Std = fit_model$tailindexStd)
  
  # Compute scale and loc for original and tail_index model estimated
  
  
  pot_measure <- PoTselect(data_measure[[col_measure]],0.1,12/timestep_measure)
  print(sum(is.na(pot_measure$pot)))
  l0_measure <- round(th_measure*length(pot_measure$pot))-1
  fit_measure <- FitGP_MLE2(pot_measure$pot, 1, N= nrow(data_measure), r11= 1, l0= l0_measure)
  fit_measure_model <- FitGP_MLE2(pot_measure$pot, 1, N= nrow(data_measure),fixedpar = fixedpar, r11= 1, l0= l0_measure)
  
  # plot ECDFs
  xrange = c(min(c(-sort(-pot_measure$pot)[l0_measure],-sort(-pot_model$pot)[l0_model])),max(c(pot_measure$pot,pot_model$pot)))
  xspan = seq(xrange[1],xrange[2],length.out=100)
  
  plot(ecdf(-sort(-pot_model$pot)),pch=1,col='black',
       xlim = xrange,
       ylim = c(1-max(th_model,th_measure),1),
       xlab = 'wind speed (m/s)',
       ylab = 'probability', 
       main = 'CDF of POT on measurements and model')
  lines(ecdf(-sort(-pot_measure$pot)),cex = 0.1,col='red')
  lines(xspan,cgp(xspan, loc = fit_model$location,scale = fit_model$scale,tail = fit_model$tailindex, t = length(pot_model$pot)/l0_model),lty='dashed',col='green')
  lines(xspan,cgp(xspan,loc = fit_measure$location, scale = fit_measure$scale, tail = fit_measure$tailindex, t = length(pot_measure$pot)/l0_measure),lty='dashed',col='blue')
  lines(xspan,cgp(xspan, loc = fit_measure_model$location, scale = fit_measure_model$scale,tail = fit_measure_model$tailindex, t = length(pot_measure$pot)/l0_measure),lty=3,col='orange')
  abline(v=c(-sort(-pot_measure$pot)[l0_measure],-sort(-pot_model$pot)[l0_model]),col=c('red','black'))
  legend('bottomright',cex = 0.6,legend=c("ecdf PoT model","ecdf PoT measurements", "cgp fit model", "cgp fit measurements","cgp fit model\n+ measurements"),
         col = c('black','red','green','blue','orange'),
         pch = c(1,1,NA,NA,NA),
         lty = c(1,1,2,2,3))
  
  print("")
  
  
  
  # mean timestep for POT select on measurements
  timestep_POT_measure <- mean(diff(pot_measure$ind)*timestep_measure/(24*365.2425)) # mean of the timediff between each peak
  
  # Return periods
  X <- 10^seq(log10(timestep_POT_measure*length(pot_measure$pot)/l0_measure), log10(10000), length.out = 30)
  
  
  print(paste('Model estimation : tail =',fit_model$tailindex,'sd',fit_model$tailindexStd,'scale =',fit_model$scale,' loc =',fit_model$location))
  # measurement original estimation
  loc = list(mean = fit_measure$location, sd = fit_measure$locationStd)
  scale = list(mean = fit_measure$scale, sd = fit_measure$logdispStd)
  tail = list(mean = fit_measure$tailindex, sd = fit_measure$tailindexStd)
  print(paste('Original estimation : tail =',tail$mean,'sd',tail$sd,'scale =',scale$mean,' loc =',loc$mean))
  
  measure_est = qgp_distrib(1-timestep_POT_measure/X,tail = tail,scale = scale,loc = loc,t=length(pot_measure$pot)/l0_measure)
  
  df_measure <- data.frame(
    return = X,
    original = measure_est$mean,
    lb = measure_est$lb,
    ub = measure_est$ub
  )
  
  
  # Estimate return values using the model-estimated tail index
  loc = list(mean = fit_measure_model$location, sd = fit_measure_model$locationStd)
  scale = list(mean = fit_measure_model$scale, sd = fit_measure_model$logdispStd)
  tail = list(mean = fit_measure_model$tailindex, sd = fit_measure_model$tailindexStd)
  print(paste('Measurements + Model estimation : tail =',tail$mean,'sd',tail$sd,'scale =',scale$mean,' loc =',loc$mean))
  
  model_est <- qgp_distrib(p = 1 - timestep_POT_measure / X, loc = loc, scale = scale, tail = tail,t=length(pot_measure$pot)/l0_measure)
  
  
  
  # Create a data frame for plotting
  df <- data.frame(
    return = X,
    original = df_measure$original,
    model_est = model_est$mean,
    lb = model_est$lb,
    ub = model_est$ub
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
    geom_line(data = df_measure, aes(y = original, color = 'all parameters estimated\nfrom measurements')) +
    geom_line(data = df_measure, aes(y = lb), lty = 'dashed', color = "red") +
    geom_line(data = df_measure, aes(y = ub), lty = 'dashed', color = "red") +
    scale_x_log10(labels = scales::comma) +
    ylim(range(c(df$original, df$model_est, df$lb, df$ub, points$speed))) +
    labs(x = 'Return Period (years)', y = 'Return Wind Speed (m/s)', title = paste0('Estimated wind speed for return period until 10000 years,\nPoT, GP, on measurements')) +
    scale_color_manual(values = c('shape estimated with model' = custom_color, 'all parameters estimated\nfrom measurements' = 'red','TRUE'='black','FALSE'='lightgrey')) +
    scale_fill_manual(values = c('Confidence Interval, 95%' = custom_color)) +
    theme(legend.position = "right")
  
  last_values <- tail(df, 1)  # Get the last row of df_measure
  print(last_values)
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
  return(list(gg = gg, df = df))
}


# wind speed estimation for POT method
# timestep in hours
# threshold : the percentage of data w want to use for the estimation
WSEst_GW <- function(data,timestep,threshold,col="F010"){
  
  pot <- PoTselect(na.approx(data[[col]]),0.1,12/timestep)
  l0 <- round(threshold*length(pot$pot))-1
  pot_fit <- FitGW_iHill(pot$pot,p=1,l0 = l0)
  
  
  timestep_selected_peaks = mean(diff(pot$ind)*timestep/(24*365.2425)) # mean of the timediff between each peak
  
  
  points <- data.frame(speed = sort(pot$pot))
  points$p <- seq(0,1-1/length(points$speed),length.out=length(points$speed))     # observed Values
  points$used <- 1-points$p <= l0/length(pot$pot)
  points$return <- timestep_selected_peaks / (1- points$p)
  
  
  print(ggplot(points,aes(x=speed,y=p,color=used))+geom_line()+geom_point()+scale_color_manual(name="",values = c('TRUE'='black','FALSE'='lightgrey')))
  
  
  tail <- list(mean=pot_fit$tailindex,sd=pot_fit$tailindexStd)
  scale <- list(mean=pot_fit$scale,sd=pot_fit$logdispStd)
  loc <- list(mean=pot_fit$location,sd=pot_fit$locationStd)
  
  X <- 10^seq(log10(timestep_selected_peaks*length(pot$pot)/l0),log10(10000),length.out = 30)
  
  
  
  distrib = qgw_distrib(1-timestep_selected_peaks/X,tail = tail,scale = scale,loc = loc,y=pot_fit$y)
  
  
  df <- data.frame(
    return = X,
    original = distrib$mean,
    lb = distrib$lb,
    ub = distrib$ub
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





WSEst_model_to_measure_GW <- function(data_model, data_measure, col_model='PF010',col_measure='F010',timestep_model=1,timestep_measure=1, th_model=NULL, th_measure=NULL, tag=NULL) {
  
  # Bootstrap on model to get uncertainties of tail index estimation
  cat("Estimation on model\n")
  flush.console()
  pot_model <- PoTselect(data_model[[col_model]],0.1,12/timestep_model)
  l0_model <- round(th_model*length(pot_model$pot))-1
  fit_model <- FitGW_iHilli(pot_model$pot,p= 1, l0= l0_model)
  
  
  
  
  # Function to constrain the tail index parameter in the estimation
  fixedpar = list(theta0 = fit_model$tailindex,theta0Std = fit_model$tailindexStd)
  
  # Compute scale and loc for original and tail_index model estimated
  
  
  pot_measure <- PoTselect(na.approx(data_measure[[col_measure]]),0.1,12/timestep_measure)
  l0_measure <- round(th_measure*length(pot_measure$pot))-1
  fit_measure <- FitGW_iHilli(pot_measure$pot,p= 1, l0= l0_measure)
  fit_measure_model <- FitGW_iHilli(pot_measure$pot,p= 1,fixedpar = fixedpar, l0= l0_measure)
  
  # plot ECDFs
  xrange = c(min(c(-sort(-pot_measure$pot)[l0_measure],-sort(-pot_model$pot)[l0_model])),max(c(pot_measure$pot,pot_model$pot)))
  xspan = seq(xrange[1],xrange[2],length.out=100)
  
  if (!is.null(tag)){
    plot(ecdf(-sort(-pot_model$pot)),pch=1,col='black',
         xlim = xrange,
         ylim = c(1-max(th_model,th_measure),1),
         xlab = 'wind speed (m/s)',
         ylab = 'probability', 
         main = 'CDF of POT on measurements and model')
    lines(ecdf(-sort(-pot_measure$pot)),cex = 0.1,col='red')
    lines(xspan,cgw(xspan, loc = fit_model$location,scale = fit_model$scale,tail = fit_model$tailindex, y= fit_model$y),lty='dashed',col='green')
    lines(xspan,cgw(xspan,loc = fit_measure$location, scale = fit_measure$scale, tail = fit_measure$tailindex,  y= fit_measure$y),lty='dashed',col='blue')
    lines(xspan,cgw(xspan, loc = fit_measure_model$location, scale = fit_measure_model$scale,tail = fit_measure_model$tailindex,  y= fit_measure_model$y),lty=3,col='orange')
    abline(v=c(-sort(-pot_measure$pot)[l0_measure],-sort(-pot_model$pot)[l0_model]),col=c('red','black'))
    legend('bottomright',cex = 0.6,legend=c("ecdf PoT model","ecdf PoT measurements", "cgw fit model", "cgw fit measurements","cgw fit model\n+ measurements"),
           col = c('black','red','green','blue','orange'),
           pch = c(1,1,NA,NA,NA),
           lty = c(1,1,2,2,3))
  }
  print("")
  
  
  
  
  # mean timestep for POT select on measurements
  timestep_POT_measure <- mean(diff(pot_measure$ind)*timestep_measure/(24*365.2425)) # mean of the timediff between each peak
  
  # Return periods
  X <- 10^seq(log10(timestep_POT_measure*length(pot_measure$pot)/l0_measure), log10(10000), length.out = 30)
  
  
  print(paste('Model estimation : tail =',fit_model$tailindex,'sd',fit_model$tailindexStd,'scale =',fit_model$scale,' loc =',fit_model$location))
  # measurement original estimation
  loc = list(mean = fit_measure$location, sd = fit_measure$locationStd)
  scale = list(mean = fit_measure$scale, sd = fit_measure$logdispStd)
  tail = list(mean = fit_measure$tailindex, sd = fit_measure$tailindexStd)
  print(paste('Original estimation : tail =',tail$mean,'sd',tail$sd,'scale =',scale$mean,' loc =',loc$mean))
  
  measure_est = qgw_distrib(1-timestep_POT_measure/X,tail = tail,scale = scale,loc = loc, y= fit_measure$y)
  
  
  
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
    model_est = model_est$mean,
    lb = model_est$lb,
    ub = model_est$ub
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
    labs(x = 'Return Period (years)', y = 'Return Wind Speed (m/s)', title = paste0('Estimated wind speed for return period until 10000 years,\nPoT + GW, on measurements')) +
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
  return(list(gg = gg, df = df,points = points))
}


#==============================================================
# GEV FIT BUT DOESN'T ACCOUNT FOR THE ESTIMATES' DEPENDENCE

WSEst_BM_2 <- function(out, bm=NULL,data=NULL,block_length=365.25, alpha = 0.05,plot_all_lines=F,col = 'F010',length.out=30,fixtailto0=F){ # we consider we discard feb 29 from datasets so no leap year
  # block_length <- as.numeric(days(block_length),'days') 
  
  if (is.null(bm)){bm <- BM_select(data,block_length,height = col)}
  
  bm <- bm$max
  
  loc = list(mean = out$original['loc'], sd = sd(out$df$loc))
  scale = list(mean = out$original['scale'], sd = sd(out$df$scale))
  tail = list(mean = out$original['tail'], sd = sd(out$df$tail))
  
  # print(bm$max)
  z <- qnorm(1 - alpha / 2)
  
  X <- 10^seq(log10(2),log10(10000),length.out = length.out)*365.25/block_length
  X <- sort(unique(c(c(50,100,1000,10000),X)))
  
  quantile_est <- qgev_distrib(exp(-1/X),tail = tail,scale = scale,loc = loc)
  print(quantile_est)
  quantile_est[quantile_est>3000] <- NA
  
  
  points <- data.frame(speed = sort(bm),return = (block_length/365.25)*1/(1-seq(1/length(bm),1-1/length(bm),length.out=length(bm))))            # observed Values
  
  
  
  
  
  if (plot_all_lines){
    
    data = data.frame(
      return = X*block_length/365.25,
      original = qgev(p=1 - 1 / X, loc = trueTail['loc'], scale = trueTail['scale'], tail = trueTail['tail']), # fit for original dataset
      boot = t(return_values)
    )    
    
    df_plot = cbind(return = X*block_length/365.25,melt(data.frame(boot=t(return_values))))
    original <- data.frame(
      return =  X*block_length/365.25, 
      original = qgev(p=1 - 1 / X, loc = trueTail['loc'], scale = trueTail['scale'], tail = trueTail['tail']) # fit for original dataset
    )
    print(df_plot)
    
    gg = ggplot(df_plot,aes(x=return)) +
      geom_line(aes(y=value,line=variable),col='gray') +
      geom_line(data = original,aes(x=return,y=original),col='red') +
      scale_x_log10()
    
  }else{
    # Create a data frame for plotting
    
    
    data <- data.frame(
      return = X*block_length/365.25,
      original = quantile_est$mean,
      lb = quantile_est$lb,
      ub = quantile_est$ub,
      sd = quantile_est$sd
    )
    
    
    # Plot using ggplot2
    gg = ggplot(data, aes(x = return)) +
      geom_ribbon(aes(ymin = lb, ymax = ub, fill = 'Confidence Interval, 95% bootstrap'), alpha = 0.2) +
      # geom_ribbon(aes(ymin = lb, ymax = ub, fill = 'Confidence Interval,  95% log(bootstrap)'), alpha = 0.2) +
      geom_line(aes(x=return,y=original,linetype = 'Original estimation'),color='red') +
      geom_point(data=points,aes(x = return,y = speed)) +
      scale_x_log10(labels= scales::comma) +
      labs(x = 'Return Period (years)', y = 'Wind Speed (m/s)', title = paste0('Estimated wind speed for return period until 10000 years, ', block_length ,'d BM (GEV)')) +
      scale_linetype_manual(name="", values = c('Original estimation' = "solid")) +
      scale_fill_manual(name = "", values = c('Confidence Interval, 95% bootstrap' = 'blue','Confidence Interval,  95% log(bootstrap)'='red')) +
      theme_minimal() +
      theme(legend.position = "right")
    
    last_values <- tail(data, 1)  # Get the last row of df_measure
    # print(last_values)
    # Annotate the last values
    gg <- gg +
      geom_text(data = last_values, aes(label = sprintf("%.2f", original), 
                                        x = max(data$return), y = original ), 
                color = 'red', size = 3, hjust = 0) +
      geom_text(data = last_values, aes(label = sprintf("%.2f", lb), 
                                        x = max(data$return), y = lb), 
                color = 'red', size = 3, hjust = 0) +
      geom_text(data = last_values, aes(label = sprintf("%.2f", ub), 
                                        x = max(data$return), y = ub), 
                color = 'red', size = 3, hjust = 0)
    
    closest_to_50 <- data %>% 
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
    
  }
  
  return(list(gg=gg, pred=data,  fit = list(tail=tail,scale=scale,loc=loc)))
}

WSEst_model_to_measure_BM_2 <- function(data_model,data_measure, col_model='PF010', col_measure='F010', Nb_model=100, Nb_measure=100, tag=NULL,length.out=30,fixtailto0=F) {
  
  if (fixtailto0){
    fit_model <- list(tail = 0,tailStd = 0)
    fit_model_pardis <- NULL
  }else{
    # Bootstrap on model to get uncertainties of tail index estimation
    cat("Bootstrap on model\n")
    flush.console()
    fit_model <- bootstrap(data_model, Nb = Nb_model, method = 'BM', column = col_model)
    fit_model_pardis <- list(tail = list(mean = fit_model$original['tail'],sd= sd(fit_model$df$tail)),
                             scale = list(mean = fit_model$original['scale'],sd= sd(fit_model$df$scale)),
                             loc = list(mean = fit_model$original['loc'],sd= sd(fit_model$df$loc)))
    fit_model <- list(tail = fit_model$original['tail'],tailStd = sd(fit_model$df$tail))
    
  }
  
  
  # Get the block maxima from measurement dataset
  bm <- BM_select(data_measure, height = col_measure)
  
  # Function to constrain the tail index parameter in the estimation
  fpar <- function(p, xpar) {
    loc <- matrix(p[1], xpar[1], xpar[2])
    scale <- matrix(p[2], xpar[1], xpar[2])
    shape <- matrix(fit_model$tail, xpar[1], xpar[2])
    list(loc = loc, scale = scale, shape = shape)
  }
  
  # Arguments for the BM_fit function
  kwargs <- list(start = c(mean(bm$max), sd(bm$max)), fpar = fpar, xpar = c(length(bm$max), 1))
  model_est <- BM_fit(bm,plot = F,kwargs = kwargs)
  
  # Quantile for confidence interval (CI)
  z <- qnorm(1 - 0.05 / 2)
  
  fit_measure <- data.frame()
  
  for (i in 1:Nb_model){
    cat('Bootstrap on Measurements:', i, '/', Nb_model, '\n')
    flush.console()
    
    
    
    # Update fpar function for each bootstrap iteration
    fpar <- function(p, xpar) {
      loc <- matrix(p[1], xpar[1], xpar[2])
      scale <- matrix(p[2], xpar[1], xpar[2])
      shape <- matrix(rnorm(n = 1,mean = fit_model$tail,sd = fit_model$tailStd), xpar[1], xpar[2])
      list(loc = loc, scale = scale, shape = shape)
    }
    
    # Arguments for the constrained bootstrap on measurements
    kwargs$fpar <- fpar
    
    # Perform bootstrap on measurement data
    fit_measure <- rbind(fit_measure,bootstrap(data_measure, Nb = Nb_measure,column = col_measure, method = 'BM', kwargs = kwargs)$df)
  }
  
  
  fit_measure <-data.frame(
    tailindex = model_est['tail'],
    tailindexStd = sd(fit_measure$tail),
    location = model_est['loc'],
    locationStd = sd(fit_measure$loc),
    scale = model_est['scale'],
    scaleStd = sd(fit_measure$scale)
  )
  
  
  # quantile estimation
  ## Return periods
  X <- 10^seq(log10(2), log10(10000), length.out = length.out)
  X <- sort(unique(c(c(50,100,1000,10000),X))) 
  
  loc = list(mean = fit_measure$location, sd = fit_measure$locationStd)
  scale = list(mean = fit_measure$scale, sd = fit_measure$scaleStd)
  tail = list(mean = fit_measure$tailindex, sd = fit_measure$tailindexStd)
  
  
  model_est = qgev_distrib(exp(-1/X),tail = tail,scale = scale,loc = loc)
  
  # Compute prediction intervals from the bootstrap measurements
  fit_measure <- WSEst_BM_2(bootstrap(data_measure, Nb = Nb_measure, column = col_measure, method = 'BM'),col = col_measure, data = data_measure,length.out=length.out)
  measure_est <-fit_measure$pred
  
  
  # Create a data frame for observed values
  points <- data.frame(speed = sort(bm$max), return = 1 / (1 - seq(1 / length(bm$max), 1 - 1 / length(bm$max), length.out = length(bm$max))))
  
  # Create a data frame for plotting
  df <- data.frame(
    return = X,
    original = measure_est$original,
    lb_o = measure_est$lb,
    ub_o = measure_est$ub,
    sd_o = measure_est$sd,
    model_est = model_est$mean,
    lb = model_est$lb,
    ub = model_est$ub,
    sd = model_est$sd,
    proba = 1/X
  )
  
  
  # Create the ggplot
  custom_color <- "blue"
  
  gg <- ggplot(df, aes(x = return)) +
    geom_line(aes(y = model_est, color = 'shape estimated with model')) +
    geom_line(aes(y = original, color = 'all parameters estimated\nfrom measurements')) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = 'Confidence Interval, 95%'), alpha = 0.2) +
    geom_line(aes(y = lb_o), lty = 'dashed', color = "red") +
    geom_line(aes(y = ub_o), lty = 'dashed', color = "red") +
    geom_point(data = points, aes(y = speed, shape = 'Measurements observed')) +
    scale_x_log10(labels = scales::comma) +
    ylim(range(c(df$original, df$model_est, df$lb, df$ub, points$speed))) +
    labs(x = 'Return Period (years)', y = 'Return Wind Speed (m/s)', title = paste0('Estimated wind speed for return period until 10000 years,\n1y BM (GEV), on measurements')) +
    scale_color_manual(values = c('shape estimated with model' = custom_color, 'all parameters estimated\nfrom measurements' = 'red')) +
    scale_fill_manual(values = c('Confidence Interval, 95%' = custom_color)) +
    theme(legend.position = "right")
  
  last_values <- tail(df, 1)  # Get the last row of df_measure
  print(last_values)
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
    geom_point(data = closest_to_50, aes(x = return, y = original), color = "red", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = model_est), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = lb), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = ub), color = "blue", size = 1) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", original), 
                                        x = return, y = original + 1), 
              color = "red", size = 3, hjust = 0) +
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
    model = fit_model_pardis,
    measure = fit_measure$fit,
    measure_model = list(tail=tail,scale=scale,loc=loc)
  )
  return(list(gg = gg, df = df,points=points, parameter_distributions = parameter_distributions))
}



#==============================================================
# GP FIT BUT DOESN'T ACCOUNT FOR THE ESTIMATES' DEPENDENCE 


# wind speed estimation for POT method
# timestep in hours
# threshold : the percentage of data w want to use for the estimation
WSEst_GP_2 <- function(data,timestep,threshold,col="F010",length.out = 30,fixtailto0=F,winter=F){
  
  ws <- na.approx(data[[col]])
  if (winter){
    ws <- na.approx(data[[col]][month(data$DateTime) %in% c(10:12,1:3)])
  }
  pot <- PoTselect_2(ws,0.3,12/timestep)
  print(length(pot$pot))
  l0 <- round(threshold*length(pot$pot))-1
  
  if (fixtailto0){fixedpar <- list(gamma0 = 0,gamma0Std = 0)}
  else{fixedpar <- NULL}
  
  pot_fit <- FitGP_MLE2(pot$pot,1, N= 0, r11= 1, fixedpar= fixedpar, l0= l0, sigma= Inf)
  
  # plot_extremes(pot,data)
  # plot(ecdf(pot$pot))
  # lines(seq(min(pot$pot),max(pot$pot),length.out=100),cgp((seq(min(pot$pot),max(pot$pot),length.out=100)-pot_fit$location)/pot_fit$scale,gamma = pot_fit$tailindex),col='red')
  
  
  timestep_selected_peaks = mean(diff(pot$ind)*timestep/(24*365.2425)) # mean of the timediff between each peak
  
  
  points <- data.frame(speed = sort(pot$pot))
  points$p <- seq(0,1-1/length(points$speed),length.out=length(points$speed))     # observed Values
  points$used <- 1-points$p <= l0/length(pot$pot)
  points$return <- timestep_selected_peaks / (1- points$p)
  
  
  # print(ggplot(points,aes(x=speed,y=p,color=used))+geom_line()+geom_point()+scale_color_manual(name="",values = c('TRUE'='black','FALSE'='lightgrey')))
  
  
  tail <- list(mean=pot_fit$tailindex,sd=pot_fit$tailindexStd)
  scale <- list(mean=pot_fit$scale,sd=pot_fit$logdispStd)
  loc <- list(mean=pot_fit$location,sd=pot_fit$locationStd)
  
  X <- 10^seq(log10(timestep_selected_peaks*length(pot$pot)/l0),log10(10000),length.out = length.out)
  X <- sort(unique(c(c(50,100,1000,10000),X)))
  
  
  distrib = qgp_distrib(1-timestep_selected_peaks/X,tail = tail,scale = scale,loc = loc,t=length(pot$pot)/l0)
  
  
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
         title = paste0('Estimated wind speed for return period until 10000 years,\nPoT: ',threshold*100,'% of the peaks, GP')) +
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





WSEst_model_to_measure_GP_2 <- function(data_model, data_measure, col_model='PF010',col_measure='F010',timestep_model=1,timestep_measure=1, th_model=NULL, th_measure=NULL, length.out = 100,peak_frac=0.3,winter=F) {
  
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
  fit_model <- FitGP_MLE2(pot_model$pot, 1, N= nrow(data_model), r11= 1, l0= l0_model)
  
  
  
  
  # Function to constrain the tail index parameter in the estimation
  fixedpar = list(gamma0 = fit_model$tailindex,gamma0Std = fit_model$tailindexStd)
  
  # Compute scale and loc for original and tail_index model estimated
  
  
  pot_measure <- PoTselect_2(ws_measure,peak_frac,12/timestep_measure)
  print(sum(is.na(pot_measure$pot)))
  l0_measure <- round(th_measure*length(pot_measure$pot))-1
  fit_measure <- FitGP_MLE2(pot_measure$pot, 1, N= nrow(data_measure), r11= 1, l0= l0_measure)
  fit_measure_model <- FitGP_MLE2(pot_measure$pot, 1, N= nrow(data_measure),fixedpar = fixedpar, r11= 1, l0= l0_measure)
  
  # # plot ECDFs
  # xrange = c(min(c(-sort(-pot_measure$pot)[l0_measure],-sort(-pot_model$pot)[l0_model])),max(c(pot_measure$pot,pot_model$pot)))
  # xspan = seq(xrange[1],xrange[2],length.out=100)
  # 
  # plot(ecdf(-sort(-pot_model$pot)),pch=1,col='black',
  #      xlim = xrange,
  #      ylim = c(1-max(th_model,th_measure),1),
  #      xlab = 'wind speed (m/s)',
  #      ylab = 'probability', 
  #      main = 'CDF of POT on measurements and model')
  # lines(ecdf(-sort(-pot_measure$pot)),cex = 0.1,col='red')
  # lines(xspan,cgp(xspan, loc = fit_model$location,scale = fit_model$scale,tail = fit_model$tailindex, t = length(pot_model$pot)/l0_model),lty='dashed',col='green')
  # lines(xspan,cgp(xspan,loc = fit_measure$location, scale = fit_measure$scale, tail = fit_measure$tailindex, t = length(pot_measure$pot)/l0_measure),lty='dashed',col='blue')
  # lines(xspan,cgp(xspan, loc = fit_measure_model$location, scale = fit_measure_model$scale,tail = fit_measure_model$tailindex, t = length(pot_measure$pot)/l0_measure),lty=3,col='orange')
  # abline(v=c(-sort(-pot_measure$pot)[l0_measure],-sort(-pot_model$pot)[l0_model]),col=c('red','black'))
  # legend('bottomright',cex = 0.6,legend=c("ecdf PoT model","ecdf PoT measurements", "cgp fit model", "cgp fit measurements","cgp fit model\n+ measurements"),
  #        col = c('black','red','green','blue','orange'),
  #        pch = c(1,1,NA,NA,NA),
  #        lty = c(1,1,2,2,3))
  # 
  # print("")
  
  
  
  # mean timestep for POT select on measurements
  timestep_POT_measure <- (length(data_measure[[col_measure]])/length(pot_measure$pot))*timestep_measure/(24*365.2425)
  
  # Return periods
  X <- 10^seq(log10(timestep_POT_measure*length(pot_measure$pot)/l0_measure), log10(10000), length.out = length.out)
  X <- sort(unique(c(c(50,100,1000,10000),X)))
  
  print(paste('Model estimation : tail =',fit_model$tailindex,'sd',fit_model$tailindexStd,'scale =',fit_model$scale,' loc =',fit_model$location))
  # measurement original estimation
  loc = list(mean = fit_measure$location, sd = fit_measure$locationStd)
  scale = list(mean = fit_measure$scale, sd = fit_measure$logdispStd)
  tail = list(mean = fit_measure$tailindex, sd = fit_measure$tailindexStd)
  print(paste('Original estimation : tail =',tail$mean,'sd',tail$sd,'scale =',scale$mean,' loc =',loc$mean))
  fit_measure = list(tail=tail,scale=scale,loc=loc)
  
  measure_est = qgp_distrib(1-timestep_POT_measure/X,tail = tail,scale = scale,loc = loc,t=length(pot_measure$pot)/l0_measure)
  
  
  
  
  # Estimate return values using the model-estimated tail index
  loc = list(mean = fit_measure_model$location, sd = fit_measure_model$locationStd)
  scale = list(mean = fit_measure_model$scale, sd = fit_measure_model$logdispStd)
  tail = list(mean = fit_measure_model$tailindex, sd = fit_measure_model$tailindexStd)
  print(paste('Measurements + Model estimation : tail =',tail$mean,'sd',tail$sd,'scale =',scale$mean,' loc =',loc$mean))
  
  model_est <- qgp_distrib(p = 1- timestep_POT_measure / X, loc = loc, scale = scale, tail = tail,t=length(pot_measure$pot)/l0_measure)
  
  
  
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
    geom_line(aes(y = original, color = 'all parameters estimated\nfrom measurements')) +
    geom_line(aes(y = lb_o), lty = 'dashed', color = "red") +
    geom_line(aes(y = ub_o), lty = 'dashed', color = "red") +
    geom_point(data = points, aes(y = speed, color=used,shape = 'Measurements observed')) +
    scale_x_log10(labels = scales::comma) +
    ylim(range(c(df$original, df$model_est, df$lb, df$ub, points$speed))) +
    labs(x = 'Return Period (years)', y = 'Return Wind Speed (m/s)', title = paste0('Estimated wind speed for return period until 10000 years,\nPoT (GP), p=',peak_frac,', th=',th_measure,', on measurements')) +
    scale_color_manual(values = c('shape estimated with model' = custom_color, 'all parameters estimated\nfrom measurements' = 'red','TRUE'='black','FALSE'='lightgrey')) +
    scale_fill_manual(values = c('Confidence Interval, 95%' = custom_color)) +
    theme(legend.position = "right")
  
  last_values <- tail(df, 1)  # Get the last row of df_measure
  print(last_values)
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
    geom_point(data = closest_to_50, aes(x = return, y = original), color = "red", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = model_est), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = lb), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = ub), color = "blue", size = 1) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", original), 
                                        x = return, y = original + 1), 
              color = "red", size = 3, hjust = 0) +
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
                 loc = list(mean = fit_model$location,sd= fit_model$locationStd)),
    measure = fit_measure,
    measure_model = list(tail=tail,scale=scale,loc=loc)
  )
  return(list(gg = gg, df = df,timestep = timestep_POT_measure,points=points, parameter_distributions = parameter_distributions))
}















