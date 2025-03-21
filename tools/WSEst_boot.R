
########################################################################################

# WIND SPEED ESTIMATION WITH BOOTSTRAP UNCERTAINTIES FOR ALL METHODS AND DISTRIBUTIONS

########################################################################################

# wind speed extrapolation using Block maxima method
# out : the output of bootstrap
WSEst_BM <- function(out, trueTail=NULL, bm=NULL,data=NULL,block_length=365.25, alpha = 0.05,plot_all_lines=F,col = 'F010',length.out=30){ # we consider we discard feb 29 from datasets so no leap year
  
  if (is.null(bm)){bm <- BM_select(data,block_length,height = col)}
  if (is.null(trueTail)){trueTail <- BM_fit(bm,block_length,plot=FALSE)}
  
  
  bm <- bm$max
  
  # print(bm$max)
  z <- qnorm(1 - alpha / 2)
  
  X <- 10^seq(log10(2),log10(10000),length.out = length.out)*365.25/block_length
  X <- sort(unique(c(c(50,100,1000,10000),X)))
  
  return_values <- matrix(NA, nrow = nrow(out$df), ncol = length(X))
  # Compute the return values
  for (i in 1:nrow(out$df)){
    return_values[i,] <- qgev(p=exp(- 1 / X), loc = out$df$loc[i], scale = out$df$scale[i], tail = out$df$tail[i])
    if (any(return_values[i,]>3000)){
      return_values[i,]<-NA
    }
  }
  
  
  points <- data.frame(speed = sort(bm),return = (block_length/365.25)*1/(1-seq(1/length(bm),1-1/length(bm),length.out=length(bm))))            # observed Values
  
  
  
  
  
  if (plot_all_lines){
    
    data = data.frame(
      return = X*block_length/365.25,
      original = qgev(p=exp( - 1 / X), loc = trueTail['loc'], scale = trueTail['scale'], tail = trueTail['tail']), # fit for original dataset
      boot = t(return_values)
    )    
    
    df_plot = cbind(return = X*block_length/365.25,melt(data.frame(boot=t(return_values))))
    original <- data.frame(
      return =  X*block_length/365.25, 
      original = qgev(p=exp( - 1 / X), loc = trueTail['loc'], scale = trueTail['scale'], tail = trueTail['tail']) # fit for original dataset
    )
    print(df_plot)
    
    gg = ggplot(df_plot,aes(x=return)) +
      geom_line(aes(y=value,line=variable),col='gray') +
      geom_line(data = original,aes(x=return,y=original),col='red') +
      scale_x_log10()
    
  }else{
    # Create a data frame for plotting
    original = qgev(p=exp(- 1 / X), loc = trueTail['loc'], scale = trueTail['scale'], tail = trueTail['tail']) # fit for original dataset
    sd = apply(log(return_values),2,sd,na.rm=TRUE)

    data <- data.frame(
      return = X*block_length/365.25,
      original = original,
      lb = original*exp(-z*sd),
      ub = original*exp(z*sd),
      sd = sd
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
  
  return(list(gg=gg, pred=data))
}



# ==============================================
# wind speed extrapolation using Block maxima method, comparing extrapolation using tail index from model simulations and from measurements. Uncertainties with bootstrap method.

WSEst_model_to_measure_BM <- function(data_model,data_measure, col_model='PF010', col_measure='F010', Nb_model=100, Nb_measure=100, tag=NULL,length.out=30, fixtailto0=F) {
  if (!fixtailto0){
  # Bootstrap on model to get uncertainties of tail index estimation
  cat("Bootstrap on model\n")
  flush.console()
  fit_model <- bootstrap(data_model, Nb = Nb_model, method = 'BM', column = col_model)
  }else{
    fit_model = list()
    fit_model$original['tail'] = 0
    fit_model$df$tail = rep(0,Nb_model)
    
  }
  # Get the block maxima from measurement dataset
  bm <- BM_select(data_measure, height = col_measure)
  
  # Function to constrain the tail index parameter in the estimation
  fpar <- function(p, xpar) {
    loc <- matrix(p[1], xpar[1], xpar[2])
    scale <- matrix(p[2], xpar[1], xpar[2])
    shape <- matrix(fit_model$original['tail'], xpar[1], xpar[2])
    list(loc = loc, scale = scale, shape = shape)
  }
  
  # Arguments for the BM_fit function
  kwargs <- list(start = c(mean(bm$max), sd(bm$max)), fpar = fpar, xpar = c(length(bm$max), 1))
  
  # Compute scale and loc for original and tail_index model estimated
  trueTail <- BM_fit(BM_select(data_measure, height = col_measure), plot = FALSE)
  modelTail <- BM_fit(BM_select(data_measure, height = col_measure), kwargs = kwargs, plot = FALSE)
  
  # Quantile for confidence interval (CI)
  z <- qnorm(1 - 0.05 / 2)
  # Return periods
  X <- 10^seq(log10(2), log10(10000), length.out = length.out)
  X <- sort(unique(c(c(50,100,1000,10000),X)))
  
  # Initialize matrix to store return values
  return_values <- matrix(NA, nrow = Nb_model * Nb_measure, ncol = length(X))
  
  for (i in 1:Nb_model) {
    cat('Bootstrap on Measurements:', i, '/', Nb_model, '\n')
    flush.console()
    
    # Update fpar function for each bootstrap iteration
    fpar <- function(p, xpar) {
      loc <- matrix(p[1], xpar[1], xpar[2])
      scale <- matrix(p[2], xpar[1], xpar[2])
      shape <- matrix(fit_model$df$tail[i], xpar[1], xpar[2])
      list(loc = loc, scale = scale, shape = shape)
    }
    
    # Arguments for the constrained bootstrap on measurements
    kwargs <- list(start = c(mean(bm$max), sd(bm$max)), fpar = fpar, xpar = c(length(bm$max), 1))
    
    # Perform bootstrap on measurement data
    fit_measure <- bootstrap(data_measure, Nb = Nb_measure,column = col_measure, method = 'BM', kwargs = kwargs)
    
    # Compute return values for each bootstrap sample
    for (j in 1:Nb_measure) {
      return_values[(i-1) * Nb_measure + j, ] <- qgev(p = exp( - 1 / X), loc = fit_measure$df$loc[j], scale = fit_measure$df$scale[j], tail = fit_measure$df$tail[j])
    }
  }
  
  # Compute standard deviation of the log of the return values
  sd <- apply(log(return_values), 2, sd,na.rm=TRUE)
  # Estimate return values using the model-estimated tail index
  model_est <- qgev(p = exp(- 1 / X), loc = modelTail['loc'], scale = modelTail['scale'], tail = modelTail['tail'])
  
  # Create a data frame for observed values
  points <- data.frame(speed = sort(bm$max), return = 1 / (1 - seq(1 / length(bm$max), 1 - 1 / length(bm$max), length.out = length(bm$max))))
  
  # Compute prediction intervals from the bootstrap measurements
  measure_est <- WSEst_BM(bootstrap(data_measure, Nb = Nb_measure, column = col_measure, method = 'BM'),col = col_measure, data = data_measure,length.out=length.out)$pred
  
  
  # Create a data frame for plotting
  df <- data.frame(
    return = X,
    original = measure_est$original,
    lb_o = measure_est$lb,
    ub_o = measure_est$ub,
    sd_o = measure_est$sd,
    model_est = model_est,
    lb = model_est * exp(-z * sd),
    ub = model_est * exp(z * sd),
    sd = sd,
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

# ==============================================
# wind speed extrapolation using POT method and GP distribution, comparing extrapolation using tail index from model simulations and from measurements. Uncertainties with bootstrap method.

WSEst_model_to_measure_GP_3 <- function(data_model, data_measure, col_model='PF010',col_measure='F010',timestep_model=1,timestep_measure=1, th_model=NULL, th_measure=NULL,Nb=10, length.out = 100,peak_frac=0.3,winter=F) {
  
  if (winter){
    ws_model <- data.frame(na.approx(data_model[month(data_model$DateTime) %in% c(10:12,1:3),c("Year",col_model)]))
    ws_measure <- data.frame(na.approx(data_measure[month(data_measure$DateTime) %in% c(10:12,1:3),c("Year",col_measure)]))
  }else{
    ws_model <- data.frame(na.approx(data_model[,c("Year",col_model)]))
    ws_measure <- data.frame(na.approx(data_measure[,c("Year",col_measure)]))
  }
  
  
  
  # Bootstrap on model to get uncertainties of tail index estimation
  cat("Estimation on model\n")
  flush.console()
  pot_model <- PoTselect_2(ws_model[[col_model]],peak_frac,12/timestep_model)
  fit_model <- bootstrap_GP(ws_model, Nb = Nb, column = col_model,timesep = 12/timestep_model,th = th_model)
  fit_model <- list(tail = fit_model$original$original$tail,tailStd = sd(fit_model$df$tail))
  
  
  
  
  
  # Compute all parameters from original measurements
  fit_measure <- bootstrap_GP(ws_measure,Nb = Nb,column = col_measure,timesep = 12/timestep_measure,th = th_measure)
  fit_measure <- data.frame(
    tail = fit_measure$original$tailindex,
    tailStd = sd(fit_measure$df$tail),
    loc = fit_measure$original$location,
    locStd = sd(fit_measure$df$loc),
    scale = fit_measure$original$scale,
    scaleStd = sd(fit_measure$df$scale)
  )
  
  # fit on measurements for tail index fixed
  pot_measure <- PoTselect_2(ws_measure[[col_measure]],peak_frac,12/timestep_measure)
  l0_measure <- round(th_measure*length(pot_measure$pot))-1
  fixedpar = list(gamma0 = fit_model$tail,gamma0Std = 0)
  fit_measure_model <- FitGP_MLE2(pot_measure$pot, 1, N= nrow(data_measure), r11= 1, l0= l0_measure,fixedpar = fixedpar)
  
  
  # uncertainties for tail index fixed
  tail_model <- rnorm(Nb,fit_model$tail, fit_model$tailStd)
  
  fit_measure_model_sd <- data.frame()
  for (i in 1:Nb){
    cat('Bootstrap on Measurements:', i, '/', Nb, '\n')
    flush.console()
    fixedpar = list(gamma0 = tail_model[i],gamma0Std = 0)
    fit_measure_model_sd <- rbind(fit_measure_model_sd,bootstrap_GP(ws_measure,Nb = Nb,column = col_measure,timesep = 12/timestep_measure,th = th_measure,fixedpar = fixedpar)$df)
  }
  
  fit_measure_model_sd = apply(fit_measure_model_sd,2,sd)
  
  fit_measure_model <- data.frame(
    tail = fit_model$tail,
    tailStd = fit_model$tailStd,
    loc = fit_measure_model$location,
    locStd = fit_measure_model_sd['loc'],
    scale = fit_measure_model$scale,
    scaleStd = fit_measure_model_sd['scale']
  )
  
  
  print(fit_measure_model)
  
  
  
  # timestep for POT select on measurements
  timestep_POT_measure <- (length(data_measure[[col_measure]])/length(pot_measure$pot))*timestep_measure/(24*365.2425)
  
  # Return periods
  X <- 10^seq(log10(timestep_POT_measure*length(pot_measure$pot)/l0_measure), log10(10000), length.out = length.out)
  X <- sort(unique(c(c(50,100,1000,10000),X)))
  
  
  
  print(paste('Model estimation : tail =',fit_model$tail,'sd',fit_model$tailStd))
  
  # measurement original estimation
  loc = list(mean = fit_measure$loc, sd = fit_measure$locStd)
  scale = list(mean = fit_measure$scale, sd = fit_measure$scaleStd)
  tail = list(mean = fit_measure$tail, sd = fit_measure$tailStd)
  print(paste('Original estimation : tail =',tail$mean,'sd',tail$sd,'scale =',scale$mean,' loc =',loc$mean))
  
  measure_est = qgp_distrib(1-timestep_POT_measure/X,tail = tail,scale = scale,loc = loc,t=length(pot_measure$pot)/l0_measure,boot=T)
  
  
  
  
  # Estimate return values using the model-estimated tail index
  loc = list(mean = fit_measure_model$loc, sd = fit_measure_model$locStd)
  scale = list(mean = fit_measure_model$scale, sd = fit_measure_model$scaleStd)
  tail = list(mean = fit_measure_model$tail, sd = fit_measure_model$tailStd)
  print(paste('Measurements + Model estimation : tail =',tail$mean,'sd',tail$sd,'scale =',scale$mean,' loc =',loc$mean))
  
  model_est <- qgp_distrib(p = 1 - timestep_POT_measure / X, loc = loc, scale = scale, tail = tail,t=length(pot_measure$pot)/l0_measure,boot=T)
  
  
  
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
  return(list(gg = gg, df = df,timestep = timestep_POT_measure,points=points, parameter_distributions = list(tail=tail,logscale=scale,loc=loc)))
}


# ==============================================
# wind speed extrapolation using POT method and GP distribution, comparing extrapolation using tail index from model simulations and from measurements. Uncertainties with bootstrap method.

WSEst_model_to_measure_GW_3 <- function(data_model, data_measure, col_model='PF010',col_measure='F010',timestep_model=1,timestep_measure=1, th_model=NULL, th_measure=NULL,Nb=10, length.out = 100,peak_frac=0.3,winter=F) {
  
  if (winter){
    ws_model <- data.frame(na.approx(data_model[month(data_model$DateTime) %in% c(10:12,1:3),c("Year",col_model)]))
    ws_measure <- data.frame(na.approx(data_measure[month(data_measure$DateTime) %in% c(10:12,1:3),c("Year",col_measure)]))
  }else{
    ws_model <- data.frame(na.approx(data_model[,c("Year",col_model)]))
    ws_measure <- data.frame(na.approx(data_measure[,c("Year",col_measure)]))
  }
  
  
  
  # Bootstrap on model to get uncertainties of tail index estimation
  cat("Estimation on model\n")
  flush.console()
  pot_model <- PoTselect_2(ws_model[[col_model]],peak_frac,12/timestep_model)
  fit_model <- bootstrap_GW(ws_model, Nb = Nb, column = col_model,timesep = 12/timestep_model,th = th_model)
  fit_model <- list(tail = fit_model$original$original$tail,tailStd = sd(fit_model$df$tail))
  
  
  
  
  # Compute scale and loc for original + uncertainties
  fit_measure <- bootstrap_GW(ws_measure,Nb = Nb,column = col_measure,timesep = 12/timestep_measure,th = th_measure)
  fit_measure <- data.frame(
    tail = fit_measure$original$tailindex,
    tailStd = sd(fit_measure$df$tail),
    loc = fit_measure$original$location,
    locStd = sd(fit_measure$df$loc),
    scale = fit_measure$original$scale,
    scaleStd = sd(fit_measure$df$scale)
  )
  
  # fit on measurements for tail index fixed
  pot_measure <- PoTselect_2(ws_measure[[col_measure]],peak_frac,12/timestep_measure)
  l0_measure <- round(th_measure*length(pot_measure$pot))-1
  fixedpar = list(theta0 = fit_model$tailindex,theta0Std = 0)
  fit_measure_model <- FitGW_iHilli(pot_measure$pot, 1, N= nrow(data_measure), r11= 1, l0= l0_measure,fixedpar = fixedpar)
  
  
  # uncertainties for tail index fixed
  tail_model <- rnorm(Nb,fit_model$tail, fit_model$tailStd)
  
  fit_measure_model_sd <- data.frame()
  for (i in 1:Nb){
    cat('Bootstrap on Measurements:', i, '/', Nb, '\n')
    flush.console()
    fixedpar = list(theta0 = tail_model[i],theta0Std = 0)
    fit_measure_model_sd <- rbind(fit_measure_model_sd,bootstrap_GW(ws_measure,Nb = Nb,column = col_measure,timesep = 12/timestep_measure,th = th_measure,fixedpar = fixedpar)$df)
  }
  
  fit_measure_model_sd = apply(fit_measure_model_sd,2,sd)
  
  fit_measure_model <- data.frame(
    tail = fit_model$tail,
    tailStd = fit_model$tailStd,
    loc = fit_measure_model$location,
    locStd = fit_measure_model_sd['loc'],
    scale = fit_measure_model$scale,
    scaleStd = fit_measure_model_sd['scale']
  )
  
  
  print(fit_measure_model)
  
  
  
  # timestep for POT select on measurements
  timestep_POT_measure <- (length(data_measure[[col_measure]])/length(pot_measure$pot))*timestep_measure/(24*365.2425)
  
  # Return periods
  X <- 10^seq(log10(timestep_POT_measure*length(pot_measure$pot)/l0_measure), log10(10000), length.out = length.out)
  X <- sort(unique(c(c(50,100,1000,10000),X)))
  
  
  
  print(paste('Model estimation : tail =',fit_model$tail,'sd',fit_model$tailStd))
  
  # measurement original estimation
  loc = list(mean = fit_measure$loc, sd = fit_measure$locStd)
  scale = list(mean = fit_measure$scale, sd = fit_measure$scaleStd)
  tail = list(mean = fit_measure$tail, sd = fit_measure$tailStd)
  print(paste('Original estimation : tail =',tail$mean,'sd',tail$sd,'scale =',scale$mean,' loc =',loc$mean))
  
  measure_est = qgw_distrib(1-timestep_POT_measure/X,tail = tail,scale = scale,loc = loc,y=log(length(pot_measure$pot)/l0_measure),boot=T)
  
  
  
  # Estimate return values using the model-estimated tail index
  loc = list(mean = fit_measure_model$loc, sd = fit_measure_model$locStd)
  scale = list(mean = fit_measure_model$scale, sd = fit_measure_model$scaleStd)
  tail = list(mean = fit_measure_model$tail, sd = fit_measure_model$tailStd)
  print(paste('Measurements + Model estimation : tail =',tail$mean,'sd',tail$sd,'scale =',scale$mean,' loc =',loc$mean))
  
  model_est <- qgw_distrib(p = 1 - timestep_POT_measure / X, loc = loc, scale = scale, tail = tail,y=log(length(pot_measure$pot)/l0_measure),boot=T)
  
  
  
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
  return(list(gg = gg, df = df,timestep = timestep_POT_measure,points=points, parameter_distributions = list(tail=tail,logscale=scale,loc=loc)))
}





###############################################################

# PARTIAL BOOTSTRAP FOR GEV DISTRIBUTION AND BM METHOD
#   Use the bootstrap to estimate the distribution of the vector of parameters 
#   (covariance matrix computed with the bootstrap and then extrapolation of the wind
#     instead of using all bootstrap estimates and then compute
#     the standard deviation of the extrapolated wind directly, better use the first one) 



WSEst_BM_3 <- function(out, bm=NULL,data=NULL,block_length=365.25, alpha = 0.05,plot_all_lines=F,col = 'F010',length.out=30){ # we consider we discard feb 29 from datasets so no leap year
  # block_length <- as.numeric(days(block_length),'days') 
  
  if (is.null(bm)){bm <- BM_select(data,block_length,height = col)}
  
  bm <- bm$max
  
  parameters = list(mean = c(out$original['tail'],out$original['loc'],out$original['scale']),cov = cov(out$df))
  
  
  
  # print(bm$max)
  z <- qnorm(1 - alpha / 2)
  
  X <- 10^seq(log10(2),log10(10000),length.out = length.out)*365.25/block_length
  X <- sort(unique(c(c(50,100,1000,10000),X)))
  
  quantile_est <- qgev_distrib_2(1-1/X,parameters)
  
  
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
  
  return(list(gg=gg, pred=data))
}

WSEst_model_to_measure_BM_3 <- function(data_model,data_measure, col_model='w10m', col_measure='F010', Nb_model=100, Nb_measure=100, tag=NULL,length.out=30) {
  
  # Bootstrap on model to get uncertainties of tail index estimation
  cat("Bootstrap on model\n")
  flush.console()
  fit_model <- bootstrap_BM(data_model, Nb = Nb_model, column = col_model)
  fit_model <- list(tail = fit_model$original['tail'],tailStd = sd(fit_model$df$tail))
  
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
  
  
  # quantile estimation
  ## Return periods
  X <- 10^seq(log10(2), log10(10000), length.out = length.out)
  X <- sort(unique(c(c(50,100,1000,10000),X))) 
  
  
  parameters = list(mean = c(model_est['tail'],model_est['loc'],model_est['scale']),cov = cov(fit_measure))
  
  
  model_est = qgev_distrib_2(1-1/X,parameters)
  
  print(model_est)
  # Compute prediction intervals from the bootstrap measurements
  measure_est <- WSEst_BM_3(bootstrap(data_measure, Nb = Nb_measure, column = col_measure, method = 'BM'),col = col_measure, data = data_measure,length.out=length.out)$pred
  
  
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
  return(list(gg = gg, df = df,points=points, parameter_distributions = parameters))
}


