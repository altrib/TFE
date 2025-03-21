

# With these functions, we can compute the GEV, GP and GW distributions (density, quantile and cumulative)
{
  gp <- function(X, loc, scale, tail) {
     X <- (X-loc)/scale
    return(ifelse(outer(tail, X, FUN = function(g, x) g != 0), 
                  outer(tail, X, FUN = function(g, x) (1 + g * x)^(-(1 + 1 / g))),
                  outer(rep(1, length(tail)), exp(-X), FUN = "*"))/scale)
  } 
  cgp <- function(X, loc, scale, tail,t=1) {
    X <- (X-loc)/scale
    return(ifelse(outer(tail, X, FUN = function(g, x) g != 0), 
                  outer(tail, X, FUN = function(g, x) 1-(1/t)*(1 + g * x)^(-(1 / g))),
                  outer(rep(1, length(tail)), 1-(1/t)*exp(-X), FUN = "*")))
  } 
  
  qgp <- function(p, loc, scale, tail,t=1) {
    return(sapply(p, function(p) ifelse(tail != 0, 
                                        loc + (scale * (t^(-tail)*(1 - p)^(-tail) - 1) / tail),
                                        loc - scale * log(t*(1 - p)))))
  }
  
  gev <- function(X, loc, scale,tail) {
    X <- (X-loc)/scale
    return(ifelse(outer(tail, X, FUN = function(g, x) g != 0), 
                  outer(tail, X, FUN = function(g, x) (1 + g * x)^(-(1 + 1 / g))*exp(-(1 + g * x)^(- 1 / g))),
                  outer(rep(1, length(tail)), exp(-X)*exp(-exp(-X)), FUN = "*"))/scale)
  }
  
  cgev <- function(X, loc, scale, tail) {
    X <- (X-loc)/scale
    return(ifelse(outer(tail, X, FUN = function(g, x) g != 0), 
                  outer(tail, X, FUN = function(g, x) exp(-(1 + g * x)^(- 1 / g))),
                  outer(rep(1, length(tail)), exp(-exp(-X)), FUN = "*")))
  }
  
  qgev <- function(p, loc, scale, tail) {
    return(sapply(p, function(p) ifelse(tail != 0, 
                                        loc + scale * ((-log(p))^(-tail)-1) / tail,
                                        loc - scale * log(-log(p)))))
  }
  
  
  gw <- function(X, loc, scale, tail,y=1) {
    X <- (X-loc)/scale
    return(ifelse(outer(tail, X, FUN = function(g, x) g != 0), 
                  outer(tail, X, FUN = function(g, x) y*exp(-y*(1 + g * x)^(1 / g))*(1 + g * x)^((1-g) / g)),
                  outer(rep(1, length(tail)), y*exp(X)*exp(-y*exp(X)), FUN = "*"))/scale)
  }
  
    
  cgw <- function(X, loc, scale, tail,y=1) {
    X <- (X-loc)/scale
    return(ifelse(outer(tail, X, FUN = function(g, x) g != 0),
                  outer(tail, X, FUN = function(g, x) 1-exp(-y*(1 + g * x)^(1 / g))),
                  outer(rep(1, length(tail)), 1-exp(-y*exp(X)), FUN = "*")))
  }
  
  qgw <- function(p, loc, scale, tail,y=1) {
    return(sapply(p, function(p) ifelse(tail != 0,
           loc + scale * (((-log(1 - p)/y)^tail - 1) / tail),
           loc + scale * log(-log(1 - p)/y))))
  }
  
}


# These functions aims to compute quantiles of the (GEV, GP and GW distributions with confidence bounds given parameters' distributions (assuming independence between parameters)

qgev_distrib = function(p,tail,scale,loc){
  if (is.numeric(tail)){tail <- list(mean=tail,sd=0)}
  if (is.numeric(scale)){scale <- list(mean=scale,sd=0)}
  if (is.numeric(loc)){loc <- list(mean=loc,sd=0)}
  
  tail_sample = rnorm(1000,tail$mean,tail$sd)
  scale_sample = rnorm(1000,scale$mean,scale$sd)
  loc_sample = rnorm(1000,loc$mean,loc$sd)
  
  quantiles_from_uncertainties = qgev(p,loc=loc_sample,scale = scale_sample,tail = tail_sample)
  
  sd = apply(log(quantiles_from_uncertainties),2,sd)
  mean = qgev(p,loc=loc$mean,scale = scale$mean,tail = tail$mean)
  lb = mean*exp(-qnorm(0.975)*sd)
  ub = mean*exp(+qnorm(0.975)*sd)
  
  return(data.frame(mean=mean,ub=ub,lb=lb,sd=sd))
}


qgp_distrib <- function(p,tail,scale,loc,t=1,boot=F){
  if (is.numeric(tail)){tail <- list(mean=tail,sd=0)}
  if (is.numeric(scale)){scale <- list(mean=scale,sd=0)}
  if (is.numeric(loc)){loc <- list(mean=loc,sd=0)}
  
  tail_sample = rnorm(1000,tail$mean,tail$sd)
  if (boot){scale_sample = rnorm(1000,log(scale$mean),scale$sd)}
  else {scale_sample = exp(rnorm(1000,log(scale$mean),scale$sd))}
  loc_sample = rnorm(1000,loc$mean,loc$sd)
  
  
  sd = apply(log(qgp(p,loc=loc_sample,scale = scale_sample,tail = tail_sample,t=t)),2,sd)
  mean = qgp(p,loc=loc$mean,scale = scale$mean,tail = tail$mean,t=t)
  lb = mean*exp(-qnorm(0.975)*sd)
  ub = mean*exp(+qnorm(0.975)*sd)
  
  return(data.frame(mean=mean,ub=ub,lb=lb,sd=sd))
}



qgw_distrib <- function(p,tail,scale,loc,y=1,boot=F){
  if (is.numeric(tail)){tail <- list(mean=tail,sd=0)}
  if (is.numeric(scale)){scale <- list(mean=scale,sd=0)}
  if (is.numeric(loc)){loc <- list(mean=loc,sd=0)}
  
  tail_sample = rnorm(1000,tail$mean,tail$sd)
  if (boot){scale_sample = rnorm(1000,log(scale$mean),scale$sd)}
  else {scale_sample = exp(rnorm(1000,log(scale$mean),scale$sd))}
  loc_sample = rnorm(1000,loc$mean,loc$sd)
  
  
  sd = apply(log(qgw(p,loc=loc_sample,scale = scale_sample,tail = tail_sample,y=y)),2,sd)
  mean = qgw(p,loc=loc$mean,scale = scale$mean,tail = tail$mean,y=y)
  lb = mean*exp(-qnorm(0.975)*sd)
  ub = mean*exp(+qnorm(0.975)*sd)
  
  return(data.frame(mean=mean,ub=ub,lb=lb,sd=sd))
}




# These functions aims to compute quantiles of the distribution with confidence bounds given parameters' distributions
# considering correlation (dependence) between parameters' estimation

qgev_distrib_2 = function(p,parameters){
  # parameter = c(tail,loc,scale)
  if (is.numeric(parameters)){parameters_sample = t(matrix(parameters,3,1000))}
  else if (parameters$cov['tail','tail']==0){
    parameters_sample = matrix(0,1000,3,dimnames = list(1:1000,c('loc','scale','tail')))
    parameters_sample[,c('loc','scale')] = rmvnorm(1000,parameters$mean[c('loc','scale')],parameters$cov[c('loc','scale'),c('loc','scale')])
  }
  else{
    parameters_sample = rmvnorm(1000,parameters$mean,parameters$cov)
  }
  colnames(parameters_sample) = names(parameters$mean)
  
  quantiles_from_uncertainties = qgev(p,loc=parameters_sample[,'loc'],scale = parameters_sample[,'scale'],tail = parameters_sample[,'tail'])
  
  
  sd = apply(log(quantiles_from_uncertainties),2,sd)
  mean = qgev(p,loc=parameters$mean['loc'],scale = parameters$mean['scale'],tail = parameters$mean['tail'])
  lb = mean*exp(-qnorm(0.975)*sd)
  ub = mean*exp(+qnorm(0.975)*sd)
  
  return(data.frame(mean=mean,ub=ub,lb=lb,sd=sd))
}


qgp_distrib_2 = function(p,parameters){
  # parameter = c(tail,loc,scale)
  if (is.numeric(parameters)){parameters_sample = t(matrix(parameters,3,100))}
  else{
    parameters_sample = rmvnorm(100,parameters$mean,parameters$cov)
  }
  
  quantiles_from_uncertainties = qgp(p,loc=parameters_sample[,2],scale = parameters_sample[,3],tail = parameters_sample[,1])
  
  sd = apply(log(quantiles_from_uncertainties),2,sd)
  mean = qgev(p,loc=parameters$mean['loc'],scale = parameters$mean['scale'],tail = parameters$mean['tail'])
  lb = mean*exp(-qnorm(0.975)*sd)
  ub = mean*exp(+qnorm(0.975)*sd)
  
  return(data.frame(mean=mean,ub=ub,lb=lb,sd=sd))
}





# plot spectrum of the data
plot_spectrum <- function(data, column = 'F010',spans=3){
  timestep = as.numeric(data$DateTime[2] - data$DateTime[1],units='hours')
  
  # print(timestep)
  # Step 3: Compute the spectrum
  spec <- spectrum(na.approx(data[[column]]),spans=spans,plot=FALSE)
  
  
  
  # Compute the frequency axis
  
  
  # Create a data frame for plotting
  spectrum_data <- data.frame(Frequency = spec$freq/timestep, Spectrum = spec$spec*timestep)
  
  
  # Plot using ggplot2
  gg = ggplot(spectrum_data, aes(x = Frequency, y = Spectrum)) +
    geom_line() +
    scale_x_continuous(limits = c(0, 0.5)) + # Plotting up to Nyquist frequency
    scale_y_log10() + 
    labs(title = "Power Spectrum of Wind Speed Data",
         x = "Frequency",
         y = "Power") +
    theme_minimal()
  return(list(spectrum_data=spectrum_data,gg=gg))
}


# Function to create a P-P plot comparing two datasets of different sizes
pp_plot_compare <- function(data1, data2,labs=NULL) {
  
  if (is.null(labs)){
    labs = list(xlab = "Empirical CDF of Dataset 1", 
                ylab = "Empirical CDF of Dataset 2", 
                main = "P-P Plot of Dataset 1 vs Dataset 2")
  }
  
  # Create empirical CDFs
  ecdf1 <- ecdf(data1)
  ecdf2 <- ecdf(data2)
  
  # Evaluate CDFs at the combined set of unique sorted values
  combined_data <- sort(unique(c(data1, data2)))
  cdf1_values <- ecdf1(combined_data)
  cdf2_values <- ecdf2(combined_data)
  
  # Plotting
  plot(cdf1_values, cdf2_values, 
       xlab = labs$xlab, 
       ylab = labs$ylab, 
       main = labs$main)
  abline(0, 1, col = "red")  # Adds a reference line
}


