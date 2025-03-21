# With this function, we can perform a bootstap EVA on the data. We randomely draw with replace a certain number of years in the available data (each years are supposed independent) and we perform the EVA on the new dataset.
# data is the dataset to which we want to perform the EVA
# Nb is the number of bootstrap we want to do
# n_draw is the number of years we want to draw from the data, default value is the number of years in the data
# column is the name of the column to which we want to do the EVA
# timesep is the time separation between the data we select by the PoT method, the separation in time is timesep*timestep (for our data it is timestep = 10min (measurements), 1h (model))
bootstrap <- function(data, Nb = 10, n_draw = NULL, column = 'F010', data_desc = NULL, method='PoT', timesep = 1 ,fraq_min=2.5, block_length=365.25,parallel=FALSE,kwargs=NULL){
  
    name_years <- unique(data$Year) # ID of years for resampling
    p <- rep(1/length(name_years), length(name_years)) # probability for resampling
    tails <- c() # to stock the results
    scale <- c() # to stock the results
    loc <- c() # to stock the results
    indices <- c()
    
    
    if (is.null(n_draw)){n_draw <- length(name_years)} # draw sample from the same length as the data
    if (is.null(data_desc)){data_desc <- paste(deparse(substitute(data)),'|',column,difftime(data$DateTime[2],data$DateTime[1],units = 'min'),'min')} # get descrition
    
    
    # Bootstrap for the PoT method
    if (method=='PoT'){ return(bootstrap_GP(data, Nb, n_draw, column = column,timesep = timesep,fraq_min=fraq_min,kwargs=kwargs))}

    # Bootstrap for block maxima method
    else if (method=='BM'){ return(bootstrap_BM(data, Nb, n_draw, column = column,kwargs=kwargs))}
      
}










bootstrap_BM <- function(data, Nb = 10, n_draw = NULL, column = 'F010', data_desc = NULL, block_length=365.25,kwargs=NULL){
  
  name_years <- unique(data$Year) # ID of years for resampling
  p <- rep(1/length(name_years), length(name_years)) # probability for resampling
  tails <- c() # to stock the results
  scale <- c() # to stock the results
  loc <- c() # to stock the results
  indices <- c()
    
    
  if (is.null(n_draw)){n_draw <- length(name_years)} # draw sample from the same length as the data
  if (is.null(data_desc)){data_desc <- paste(deparse(substitute(data)),'|',column,difftime(data$DateTime[2],data$DateTime[1],units = 'min'),'min')} # get descrition
  

  bm <- BM_select(data,block_length,height = column)
  # bootstrap
  for (bi in 1:Nb) {
    cat("\rRun",bi,'/',Nb)
    flush.console()
    
    # bootstrap resampling
    if (is.Date(name_years)){name_years <- year(name_years)}
    index <- sample(name_years, size = n_draw, replace = TRUE, prob = p)
    indices <- rbind(indices,index)
    temp <- data.frame()
    for (d in index){
      temp <- rbind(temp,bm[bm$year == d, ])
    }
    
    # GEV distribution fitting on the reseampled dataset
    boot_fit <- BM_fit(temp,block_length,plot = F,kwargs = kwargs)
    
    
    # result treatment
    tails <- rbind(tails,boot_fit['tail'])
    scale <- rbind(scale,boot_fit['scale'])
    loc <- rbind(loc,boot_fit['loc'])
  }
  df <- data.frame(tail = tails,loc = loc,scale = scale)
  
  original <- BM_fit(bm,block_length,plot = F,kwargs = kwargs)

  cat('\n')
  
  
  # Description of the result
  desc = paste("Result for",Nb,"estimations, n_draw=",n_draw,"years,\nmethod: BM 1y")

  desc = paste(desc,"\ndata :",data_desc)
  
  desc=list(fulltext=desc,Nb=Nb,n_draw=n_draw,method="BM",timestep=difftime(data$DateTime[2],data$DateTime[1],units = 'min'))
  
  return(list(df =df,original= original,desc =  desc,years = indices))
}


bootstrap_GP <- function(data, Nb = 10, n_draw = NULL, column = 'F010', timesep = 1 ,th=NULL,fraq_min=2.5,fixedpar=NULL){
  
  name_years <- unique(data$Year) # ID of years for resampling
  p <- rep(1/length(name_years), length(name_years)) # probability for resampling
  tail <- c() # to stock the results
  scale <- c() # to stock the results
  loc <- c() # to stock the results
  indices <- c()
  sample_fraction = 10^(seq(-fraq_min, 0, length.out=30))
  
  
  if (is.null(n_draw)){n_draw <- length(name_years)} # draw sample from the same length as the data


  error = 0
  original_pot <- data.frame(PoTselect_2(na.approx(data[[column]]),0.3,timesep))
  first_occurrences <- which(!duplicated(data$Year))
  pot$Year <- sapply(original_pot$ind,function (idx){
    data$Year[max(first_occurrences[first_occurrences <= idx])]
  })
  pot_by_year <- split(original_pot, pot$Year)
  # data <- data.frame(na.approx(data[month(data$DateTime) %in% c(10:12,1:3),c('Year',column)]))
  # data_by_year <- split(data, data$Year)

  
  for (bi in 1:Nb) {
    cat("\rRun",bi,'/',Nb)
    flush.console()
    
    # bootstrap resampling 
    index <- sample(name_years, size = n_draw, replace = TRUE, prob = p)
    indices <- rbind(indices,index)
    
    # new_data <- do.call(rbind, lapply(index, function(d) data_by_year[[as.character(d)]]))
    # pot <- PoTselect_2(new_data[[column]],p = 0.3,separation = timesep)
    # rm(new_data)
    
    pot <- do.call(rbind, lapply(index, function(d) pot_by_year[[as.character(d)]]))

    
    if(!is.null(th)){l0 <- round(th*length(pot$pot))-1}
    else{l0 <- round(sample_fraction * length(pot$pot)) - 1}
    # GP distribution fitting
    try({
      out <- FitGP_MLE2(pot$pot, 1, N = 0, r11 = 1, fixedpar = fixedpar, l0 = l0, sigma = Inf, metadata = NULL)
    }, silent = TRUE)
    
    if (!exists("out")) {
      stop("Failed to fit : try to decremente fraq_min.")
    }
    
    tail <- rbind(tail,cbind(out$tailindex,rep(bi,length(out$l))))
    loc <- c(loc,out$location)
    scale <- c(scale,out$scale)
  }
  
  df <- cbind(out$l/out$N,loc,scale,tail)
  df <- data.frame(df)
  colnames(df) <- c("p","loc","scale","tail","boot")
  
  X <- PoTselect_2(na.approx(data[[column]]),0.3,timesep)
  if(!is.null(th)){l0 <- round(th*length(X$pot))-1}
  else{l0 <- round(sample_fraction * length(X$pot)) - 1}
  original <- FitGP_MLE2(X$pot, 1, N= 0, r11= 1, fixedpar= NULL, l0= l0, sigma= Inf, metadata=data.frame(label=column))
  original$original <- data.frame(p = original$l/original$N,tail = original$tailindex)
  
  print(error)

  cat('\n')
  
  return(list(df =df,original= original))
}


bootstrap_GW <- function(data, Nb = 10, n_draw = NULL, column = 'F010', timesep = 1 ,th=NULL,fraq_min=2.5,fixedpar=NULL){
  
  name_years <- unique(data$Year) # ID of years for resampling
  p <- rep(1/length(name_years), length(name_years)) # probability for resampling
  tail <- c() # to stock the results
  scale <- c() # to stock the results
  loc <- c() # to stock the results
  indices <- c()
  sample_fraction = 10^(seq(-fraq_min, 0, length.out=30))
  
  
  if (is.null(n_draw)){n_draw <- length(name_years)} # draw sample from the same length as the data
  
  
  error = 0
  original_pot <- data.frame(PoTselect_2(na.approx(data[[column]]),0.3,timesep))
  first_occurrences <- which(!duplicated(data$Year))
  pot$Year <- sapply(original_pot$ind,function (idx){
    data$Year[max(first_occurrences[first_occurrences <= idx])]
  })
  pot_by_year <- split(original_pot, pot$Year)
  # data <- data.frame(na.approx(data[month(data$DateTime) %in% c(10:12,1:3),c('Year',column)]))
  # data_by_year <- split(data, data$Year)
  
  
  for (bi in 1:Nb) {
    cat("\rRun",bi,'/',Nb)
    flush.console()
    
    # bootstrap resampling 
    index <- sample(name_years, size = n_draw, replace = TRUE, prob = p)
    indices <- rbind(indices,index)
    
    # new_data <- do.call(rbind, lapply(index, function(d) data_by_year[[as.character(d)]]))
    # pot <- PoTselect_2(new_data[[column]],p = 0.3,separation = timesep)
    # rm(new_data)
    
    pot <- do.call(rbind, lapply(index, function(d) pot_by_year[[as.character(d)]]))
    
    
    if(!is.null(th)){l0 <- round(th*length(pot$pot))-1}
    else{l0 <- round(sample_fraction * length(pot$pot)) - 1}
    # GP distribution fitting
    try({
      out <- FitGW_iHilli(pot$pot, 1, N = 0, r11 = 1, fixedpar = fixedpar, l0 = l0, sigma = Inf, metadata = NULL)
    }, silent = TRUE)
    
    if (!exists("out")) {
      stop("Failed to fit : try to decremente fraq_min.")
    }
    
    tail <- rbind(tail,cbind(out$tailindex,rep(bi,length(out$l))))
    loc <- c(loc,out$location)
    scale <- c(scale,out$scale)
  }
  
  df <- cbind(out$l/out$N,loc,scale,tail)
  df <- data.frame(df)
  colnames(df) <- c("p","loc","scale","tail","boot")
  
  X <- PoTselect_2(na.approx(data[[column]]),0.3,timesep)
  if(!is.null(th)){l0 <- round(th*length(X$pot))-1}
  else{l0 <- round(sample_fraction * length(X$pot)) - 1}
  original <- FitGW_iHilli(X$pot, 1, N= 0, r11= 1, fixedpar= NULL, l0= l0, sigma= Inf, metadata=data.frame(label=column))
  original$original <- data.frame(p = original$l/original$N,tail = original$tailindex)
  
  print(error)
  
  cat('\n')
  
  return(list(df =df,original= original))
}
