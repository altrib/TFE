library(EVTools)
library(ggplot2)
library(reshape2)
library(lubridate)
library(boot)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(evd)
library(stats)
library(zoo)
library(ggnewscale)
library(gganimate)
library(rlang)
library(plotly)
library(R.matlab)
library(akima)
library(ggstance)
library(spam)
library(ComplexHeatmap)
library(circlize)

# 
# setwd(utils::getSrcDirectory(function(){})[1])
source('../tools/BM.R')
source('../tools/distributions.R')
source('../tools/fgev.flex.R')
source('../tools/tailindexplots.R')
source('../tools/WSEst_boot.R')
source('../tools/WSEst_Analytic.R')
source('../tools/WSEst_DEPRECIATED.R')
source('../tools/bootstrap.R')
source("../tools/heatmap_plot.R")

# With this function, we can plot the extreme values selected by the PoT and BM method
plot_extremes <- function(pot,data,ylim=NULL,xlim = NULL,col='F010',m='POT',winter=NULL){
  
  if (is.null(winter)){winter = rep(TRUE,nrow(data))}
  else{
    index_winter = rle(winter)
    index_winter$cumsum <- cumsum(index_winter$lengths)
    start = 1
    for (i in seq_along(index_winter$lengths)){
      if (!index_winter$values[i]){
        pot$ind[start:length(pot$ind)] = pot$ind[start:length(pot$ind)] + index_winter$lengths[i]
      }
      else{
        start = min(which(pot$ind > index_winter$cumsum[i]))
        if (start==Inf){break}
      }
    }
  }
  
    gg = ggplot(data,aes(x=DateTime,y=.data[[col]],colour=winter)) +
      geom_line()+
      geom_point(data = data[pot$ind, ], aes(colour = "Selected"), size = 1,pch=19) +
      labs(x = "Date", y = "Wind speed (m/s)", title = paste("Selected",m)) +
      theme_minimal()+
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
    
  if (m=='POT'){
    gg = gg + geom_hline(yintercept = min(pot$pot), linetype = "dashed", color = "red") + # Horizontal line
    scale_colour_manual(values = c("Selected" = "red","TRUE"='darkblue','FALSE'='lightgrey')) # Specify red color for points
  }
  else if (m=='BM'){
    gg = gg + geom_vline(xintercept = data$DateTime[month(data$DateTime) == 7 & day(data$DateTime) == 1 & hour(data$DateTime) == 0 & minute(data$DateTime) == 0],lty='dashed',col='grey')
  }
  if (!is.null(xlim)) {
    gg = gg + xlim(xlim)
  }
  if (!is.null(ylim)){
    gg = gg + ylim(ylim)
  }else{
    gg = gg + ylim(c(0.50*min(data[[col]][pot$ind]),max(data[[col]])))
  }
  return(gg)
}






plot_par<- function(boot) {
  tail <- c(boot$original['tail'],boot$original['tail']-z*sd(boot$df$tail),boot$original['tail']+z*sd(boot$df$tail))
  loc <- c(boot$original['loc'],boot$original['loc']-z*sd(boot$df$loc),boot$original['loc']+z*sd(boot$df$loc))
  scale <- c(boot$original['scale'],boot$original['scale']-z*sd(boot$df$scale),boot$original['scale']+z*sd(boot$df$scale))
  
  plot(1,tail[1],ylim=tail[2:3])
  arrows(1,tail[2],1,tail[3],code=3,angle=90)
  
  plot(1,scale[1],ylim=scale[2:3])
  arrows(1,scale[2],1,scale[3],code=3,angle=90)
  
  plot(1,loc[1],ylim=loc[2:3])
  arrows(1,loc[2],1,loc[3],code=3,angle=90)
  return(list(tail=tail,loc=loc,scale=scale))
}





# Fonction pour une transformation pseudo-logarithmique avec un seuil ajusté
pseudo_log <- function(x, base = 10, threshold = 0.1) {
  sign(x) * log10(abs(x) / threshold + 1)
}


mean_residual_life <- function(data, thresholds) {
  mrl_values <- sapply(thresholds, function(threshold) {
    exceedances <- data[data > threshold] - threshold
    if (length(exceedances) > 0) {
      mean(exceedances)
    } else {
      NA
    }
  })
  return(mrl_values)
}



PoTselect_2 <- function(data, p=NULL, separation=1, threshold=NULL) {
  n <- length(data)
  sdata <- -sort(-data)
  
  if(!is.null(p)){
  threshold <- sdata[round(n * p)]
  # print(threshold)
  } 
  above_threshold <- data > threshold
  
  # Find the start and end indices of each excursion
  excursion_starts <- which(diff(c(0, above_threshold)) == 1)
  excursion_ends <- which(diff(c(above_threshold, 0)) == -1)

  # Ensure they match up
  if (length(excursion_starts) != length(excursion_ends)) {
    stop("Mismatch in excursion start and end points")
  }
  
  
  test_time = now()
  # Find the start and end indices of each excursion
  ind <- sapply(1:length(excursion_starts), function(i) {
    excursion_starts[i]+which.max(data[excursion_starts[i]:excursion_ends[i]])-1
  })

  # Initial potential excursions above threshold
  pot <- data[ind]
  count <- 0
  
  # Remove excursions that are too close
  repeat {
    di <- diff(ind)
    if (min(di) >= separation) {
      break
    }
    
    ldi <- length(di)
    imins <- which((di < separation) & (diff(pot) > 0))
    iplus <- which((di < separation) & (diff(pot) <= 0))
    
    pot[imins] <- 0
    pot[iplus + 1] <- 0
    
    id <- pot > 0
    ind <- ind[id]
    pot <- data[ind]
    
    count <- count + 1
  }
  
  potdata <- list("ind" = ind, "pot" = pot, "p" = p)
  
  return(potdata)
}



find_max_window <- function(data, pot_ind, window_size) {
  best_window <- NULL
  counts <- rep(0, length(data))
  for (i in seq(1,length(data),max(1,round(window_size/10)))) {
    cat("\r",i,'/',length(data))
    flush.console()
    start <- max(1, i - window_size + 1)
    end <- min(length(data), i)
    counts[i] <- sum(pot_ind>start & pot_ind < end)
    
    max_count <- max(counts)
    best_window <- c(which.max(counts),which.max(counts)+window_size)
  }
  
  
  return(list(best_window = best_window, max_count = max_count))
}




