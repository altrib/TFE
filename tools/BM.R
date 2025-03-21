# This function fits the GEV (Generalized Extreme Value) distribution 
# to the block maxima of a given data series and plots the ECDF and histogram.
BM_fit <- function(bm, block_length=365.25, plot=T, kwargs=NULL) {
  # Create a matrix of block maxima
  temp <- matrix(bm$max)
  
  # If no additional parameters are provided, define default parameters
  if (is.null(kwargs)){
    start <- c(mean(temp), sd(temp), 0)  # Starting values for loc, scale, and shape
    xpar <- c(length(temp), 1)
    
    # Function to structure parameters for the GEV fitting process
    fpar <- function(p, xpar) {
      loc <- matrix(p[1], xpar[1], xpar[2])
      scale <- matrix(p[2], xpar[1], xpar[2])
      shape <- matrix(p[3], xpar[1], xpar[2])
      list(loc = loc, scale = scale, shape = shape)
    }
  } else {
    start <- kwargs$start
    xpar <- kwargs$xpar
    fpar <- kwargs$fpar
  }
  
  # Fit the GEV distribution using maximum likelihood estimation
  out <- fgev.flex(temp, start, fpar, xpar)
  
  # Ensure the fitting process returns the correct number of parameters
  if (length(out$estimate) != 3) {
    out$estimate <- unlist(fpar(out$estimate, c(1, 1)))
  }
  
  # Extract the tail, location, and scale parameters from the fit
  trueTail = c(out$estimate[3], out$estimate[1], out$estimate[2])
  names(trueTail) = c('tail', 'loc', 'scale')
  
  # Generate a range of values for plotting
  xspan <- seq(min(temp), max(temp), length.out = 100)
  
  # If plotting is enabled, plot the histogram and ECDF
  if (plot) {
    # Plot histogram with GEV density curve
    hist(temp, breaks=10, freq=F, main=paste0('extreme wind speed distribution, 10m,\n', block_length, 'd BM'), xlab='wind speed (m/s)')
    lines(xspan, gev(xspan, loc=trueTail['loc'], scale=trueTail['scale'], tail=trueTail['tail']), col='red')
    
    # Plot ECDF with GEV cumulative distribution function
    plot(ecdf(temp), main=paste0('extreme wind speed cumulative distribution, 10m,\n', block_length, 'd BM'), xlab='wind speed (m/s)')
    lines(xspan, cgev(xspan, loc=trueTail['loc'], scale=trueTail['scale'], tail=trueTail['tail']), col='red')
  }
  
  # Return the fitted GEV parameters
  return(trueTail)
}

# This function fits the GEV distribution with covariance computation and plots it.
BM_fit_cov <- function(bm, block_length=365.25, plot=T, kwargs=NULL, cov_type='num') {
  temp <- matrix(bm$max)
  
  if (is.null(kwargs)){
    start <- c(mean(temp), sd(temp), 0)  # Default starting values for fitting
    xpar <- c(length(temp), 1)
    fpar <- function(p, xpar) {
      loc <- matrix(p[1], xpar[1], xpar[2])
      scale <- matrix(p[2], xpar[1], xpar[2])
      shape <- matrix(p[3], xpar[1], xpar[2])
      list(loc = loc, scale = scale, shape = shape)
    }
  } else {
    start <- kwargs$start
    xpar <- kwargs$xpar
    fpar <- kwargs$fpar
  }
  
  # Fit the GEV distribution
  out <- fgev.flex(temp, start, fpar, xpar)
  
  # Ensure parameters are structured properly
  if (length(out$estimate) != 3) {
    out$estimate <- unlist(fpar(out$estimate, c(1, 1)))
  }
  
  # Extract the fitted location, scale, and tail parameters
  trueTail = c(out$estimate[1], out$estimate[2], out$estimate[3])
  names(trueTail) = c('loc', 'scale', 'tail')
  
  # If analytic covariance calculation is requested, use the gevinfom function
  if (cov_type == 'analytic') {
    out$cov = gevinfom(trueTail['tail'], trueTail['scale'], length(temp))$cov
  }
  
  # If extra parameters (kwargs) are provided, adjust the covariance matrix
  if (!is.null(kwargs)) {
    cov = matrix(0, 3, 3)
    colnames(cov) <- c('loc', 'scale', 'tail')
    rownames(cov) <- c('loc', 'scale', 'tail')
    cov[c('loc', 'scale'), c('loc', 'scale')] = out$cov[c('loc', 'scale'), c('loc', 'scale')]
    cov['tail', 'tail'] = kwargs$tailStd
  } else {
    cov = out$cov
    colnames(cov) <- c('loc', 'scale', 'tail')
    rownames(cov) <- c('loc', 'scale', 'tail')
  }
  
  # Print the covariance matrix
  print(cov)
  
  # Plot if enabled
  if (plot) {
    xspan <- seq(min(temp), max(temp), length.out = 100)
    hist(temp, breaks=10, freq=F, main=paste0('extreme wind speed distribution, 10m,\n', block_length, 'd BM'), xlab='wind speed (m/s)')
    lines(xspan, gev(xspan, loc=trueTail['loc'], scale=trueTail['scale'], tail=trueTail['tail']), col='red')
    
    # Plot ECDF
    plot(ecdf(temp), main=paste0('extreme wind speed cumulative distribution, 10m,\n', block_length, 'd BM'), xlab='wind speed (m/s)')
    lines(xspan, cgev(xspan, loc=trueTail['loc'], scale=trueTail['scale'], tail=trueTail['tail']), col='red')
  }
  
  # Return the fitted parameters and covariance matrix
  return(list(trueTail=trueTail, cov=cov))
}

# Function to select block maxima from a data series based on block length
BM_select <- function(data, block_length=365.25, height='F010') {
  if (is.Date(data$Year)) {
    data$Year <- year(data$Year)
  }
  
  # Select block maxima based on yearly blocks
  if (block_length == 365.25) {
    bm <- data %>%
      mutate(row_index = row_number()) %>%
      group_by(year = Year) %>%
      summarise(
        max = max(!!sym(height), na.rm = TRUE),
        ind = row_index[which.max(!!sym(height))]
      ) %>%
      ungroup()
  } else {
    # Remove leap day and calculate block maxima
    bm <- data %>%
      filter(!(month(DateTime) == 2 & day(DateTime) == 29)) %>%
      mutate(
        block = (as.numeric(as.Date(DateTime, format="%Y-%m-%d") - ymd(Year)) %/% block_length) + 1
      ) %>%
      group_by(year=year(Year), block) %>%
      summarise(max = max(!!sym(height), na.rm = TRUE), .groups='drop') %>%
      ungroup()
  }
  
  # Return the block maxima
  return(bm)
}

# Function to plot tail index estimation across multiple heights with confidence intervals
BM_TIplot <- function(data, Nbs=c(500), heights=NULL, title=NULL, parameter='tail') {
  if (is.null(heights)) {
    heights = colnames(data)[grepl("^[Ffw]", names(data))]
  }
  
  datasets <- list()
  data_desc <- paste(deparse(substitute(data)), '|', difftime(data$DateTime[2], data$DateTime[1], units='min'), 'min')
  
  for (Nb in Nbs) {
    outs <- list()
    for (height in heights) {
      print(height)
      tmp <- bootstrap(data, Nb=Nb, column=height, method='BM', data_desc=data_desc)
      outs[[as.character(height)]] <- tmp$df[[parameter]]
      if (is.null(title)) {
        title <- tmp$desc$fulltext
      }
    }
    datasets[[as.character(Nb)]] <- data.frame(outs)
  }
  
  colors <- if (length(Nbs) > 2) {
    brewer.pal(n=length(Nbs), name="Set1")
  } else {c('red', 'blue')}
  
  offset <- 0.15
  alpha <- 0.05
  z <- qnorm(1 - alpha / 2)
  
  tailindex <- c()
  for (height in heights) {
    tailindex <- c(tailindex, BM_fit(BM_select(data, height=height), plot=F)[parameter])
  }
  
  # Set up plot with extra space for the legend
  par(mar=c(5, 4, 4, 15) + 0.1)
  
  # Initialize plot
  plot(x=1:length(heights), y=tailindex, col="blue", pch=19, cex=1,
       xlab="heights", ylim=c(min(unlist(lapply(datasets, function(x) tailindex - z * apply(x, 2, sd)))), 
                              max(unlist(lapply(datasets, function(x) tailindex + z * apply(x, 2, sd))))),
       xlim=c(1-(length(Nbs)-1)*offset, length(heights)),
       xaxt="n", ylab=parameter)
  axis(1, at=1:length(heights), labels=heights, tick=TRUE)
  
  for (i in seq_along(Nbs)) {
    Nb <- Nbs[i]
    dataset <- datasets[[as.character(Nb)]]
    color <- colors[i]
    
    lower_bounds <- tailindex - z * apply(dataset, 2, sd)
    upper_bounds <- tailindex + z * apply(dataset, 2, sd)
    
    arrows(x0=1:ncol(dataset)-offset*(i-1), y0=lower_bounds, x1=1:ncol(dataset)-offset*(i-1), 
           y1=upper_bounds, code=3, angle=90, length=0.1, col=color)
  }
  
  legend("topright", inset=c(-0.8, 0), legend=c("Value without bootstrap", paste("Nb =", Nbs)), 
         col=c("blue", colors), pch=c(19, rep(1, length(Nbs))), bty='n', xpd=TRUE)
  
  title(title)
  return(datasets)
}
