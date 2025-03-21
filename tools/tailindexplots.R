# Stability plot from EVTools but available for one or more estimations and for all parameters (tail index, scale and location)
# list_es : list of the estimations with PoT selection on GW or GP (take only list, even with only one estimation)
# overlap : if TRUE, make the stability plot from all estimation on the sam plot, if FALSE, make separated figures
# ncol : number of column if there are separated figures (overlap = FALSE)
# fix_y : if TRUE fix the y axis to the biggest span. If numeric, ylim = fix_y.
# selected_threshold : print a line to the value of this parameter on the plot and return the estimation for this threshold (should be a sample fraction). If not numeric, select a threshold where the difference between two estimation is minimum.
# parameter : Plot the stability of "tail" index, "scale", "location" or "all".
tailindexplot_2 <- function(list_es= NULL, overlap = FALSE,ncol = 3, fix_y = FALSE,selected_threshold = NULL,parameter='tail') {
  data = data.frame()
  group = c()
  selected <- NULL
  for (i in 1:length(list_es)){
    es <- list_es[[i]]
    metadata <- es$metadata
    label <- metadata$label
    if (length(label)==0) {label <- i}
    
    pconf <- 0.95
    qn <- abs(qnorm((1+pconf)/2)) # half width of normal confidence interval
    
    if (parameter=='tail'){
      ub <- es$tailindex + es$tailindexStd*qn
      lb <- es$tailindex - es$tailindexStd*qn
      data <- rbind(data,cbind(list_es[[i]]$l/list_es[[i]]$N,es$tailindex,lb,ub))
    }else if (parameter == 'scale'){
      ub <- es$scale*exp(es$logdispStd*qn)
      lb <- es$scale*exp(-es$logdispStd*qn)
      data <- rbind(data,cbind(list_es[[i]]$l/list_es[[i]]$N,es$scale,lb,ub))
    }else if (parameter == 'loc'){
      ub <- es$location + es$locationStd*qn
      lb <- es$location - es$locationStd*qn
      data <- rbind(data,cbind(list_es[[i]]$l/list_es[[i]]$N,es$location,lb,ub))
    }else if (parameter == 'logdisp'){
      ub <- es$logdisp + es$logdispStd*qn
      lb <- es$logdisp -es$logdispStd*qn
      data <- rbind(data,cbind(list_es[[i]]$l/list_es[[i]]$N,es$logdisp,lb,ub))
    }else if (parameter == 'testlogdisp'){
      ub <- exp(es$logdisp + es$logdispStd*qn)*es$location
      lb <- exp(es$logdisp -es$logdispStd*qn)*es$location
      data <- rbind(data,cbind(list_es[[i]]$l/list_es[[i]]$N,exp(es$logdisp)*es$location,lb,ub))
    }else if (parameter=='all'){
      tail = tailindexplot_2(list_es,overlap,ncol,fix_y,selected_threshold,parameter='tail')
      scale = tailindexplot_2(list_es,overlap,ncol,fix_y,selected_threshold,parameter='scale')
      loc = tailindexplot_2(list_es,overlap,ncol,fix_y,selected_threshold,parameter='loc')
      return(list(tail = tail,scale=scale, loc = loc))
    }
    group = c(group,rep(label,length(es$tailindex)))
  }
  # print(group)
  data$group = factor(group,levels=unique(group))
  colnames(data)<- c('p','val','lb','ub','group')
  # print(data)
  
  
  # axis labels 
  ylab <- paste(parameter)
  xlab <- paste("sample fraction for quantile estimation")
  
  
  options(repr.plot.width=20, repr.plot.height=12)
  
  # plot
  gg = ggplot(data, aes(x = p, y = val, color = group)) +
    geom_point() +
    geom_line() +
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha=0.02) +
    labs(x = xlab, y = ylab, color = "group") +
    scale_x_log10()
  
  if (!is.null(selected_threshold)){
    if (is.numeric(selected_threshold)){
      index = which.min(abs(data$p - selected_threshold))
      gg = gg + geom_vline(xintercept = selected_threshold)
      selected = c(tail = data$val[index],sd = es$tailindexStd[index], lb = data$lb[index], ub = data$ub[index])
    }
    else if (selected_threshold){
      index = c(which(abs(diff(es$tailindex))<10^(-4))+1,which.min(abs(diff(es$tailindex)))+1)
      gg = gg + geom_vline(xintercept = data$p[index])
      index = max(index)
      selected <- c(tail = data$val[index],sd = es$tailindexStd[index], lb = data$lb[index], ub = data$ub[index])
      print(paste0("Selected : p = ",data$p[index],selected))
    }
  }
  
  if (!overlap){gg = gg + facet_wrap(~group,ncol = ncol, scales = "free_y")+ theme(legend.position = "none")}
  if(is.numeric(fix_y)){
    gg = gg + ylim(fix_y)
  }
  else if (fix_y){
    scale_view_offset = (max(data$ub) - min(data$lb))/100
    gg = gg + ylim(min(data$lb)-scale_view_offset,max(data$ub)+scale_view_offset)
  }
  

  return(list(gg=gg,df=data,selected = selected))
  
} # klaar



#we define here a new tail index plot function to plot the result with the bootstrap uncertainties
# df : data frame output of the bootstrap function
# original : original fit from the point estimate
tailindexplot_3 <- function(df,original,pconf=0.95,desc = NULL,plot_outliers=T,fix_y=NULL){
  z <- abs(qnorm((1-pconf)/2))
  
  # Calculate quartiles for each value of p
  intervals_df <- df %>%
    group_by(p) %>%
    summarise(
      Q1 = quantile(tail, 0.25),
      Q2 = median(tail),
      Q3 = quantile(tail, 0.75),
      mean = mean(tail),
      std = sd(tail),
    ) %>%
    ungroup()
  std = -cummax(-intervals_df$std)
  
  
  intervals_df[c('lb','ub')] <- cbind(original$tail-z*std,original$tail+z*std)
  print(intervals_df)
  
  # outliers
  out_idx = df$tail<intervals_df$lb | df$tail>intervals_df$ub
  out_df = data.frame(p =df$p[out_idx],tail =df$tail[out_idx])
  
  
  # plot
  if (is.null(desc)){title = paste("Result for", n_distinct(df$boot) ,"samples")}
  else{title = desc}
  
  
  gg = ggplot(intervals_df) +
    geom_line(aes(x = p, y = mean, linetype = "Mean"), color = "black") +
    geom_line(aes(x = p, y = Q2, linetype = "Median"), color = "black") +
    # geom_ribbon(aes(x = p, ymin = Q1, ymax = Q3, fill = "1st and 3rd Quartiles"), alpha = 0.3) +
    geom_ribbon(aes(x = p, ymin = lb, ymax = ub, fill = "pconf interval"), alpha = 0.3) +
    geom_line(data=original,aes(x=p,y=tail,color="Original"))+
    geom_point(data=original,aes(x=p,y=tail,color="Original"))+
    scale_x_log10() +
    labs(
      y = "tail index",
      x = "sample fraction for quantile estimation",
      title = title,
    ) +
    scale_linetype_manual(
      name = "Legend",
      values = c("Mean" = "dashed", "Median" = "solid"),
      labels = c("Mean", "Median")
    ) +
    scale_fill_manual(
      name = "Legend",
      values = c("pconf interval"  = "blue"),
      labels = c( paste0(100*pconf,"% Interval"))
    ) +
    scale_color_manual(name = "Legend",values=c("Original" = "red"), labels=c("Original"="Original dataset"))
  if (plot_outliers){
    gg =gg + 
      geom_point(data = out_df, aes(x = p, y = tail, color = "Outliers"), size = 2) +
      scale_color_manual(
        name = "Legend",
        values = c("Outliers" = "black"),
        labels = c("Outliers")
      )}
  if(is.numeric(fix_y)){
    gg = gg + ylim(fix_y)+ theme(legend.position = "none")
  }
  
  return(gg)
}


tailindexplot_4 <- function(df,original,desc = NULL,fix_y=NULL){
# from bootstrap plot all est

  # Calculate quartiles for each value of p
  intervals_df <- df %>%
    group_by(p) %>%
    summarise(
      Q1 = quantile(tail, 0.25),
      Q2 = median(tail),
      Q3 = quantile(tail, 0.75),
      mean = mean(tail),
      std = sd(tail),
    ) %>%
    ungroup()
  
  
  # outliers
  out_idx = df$tail<intervals_df$lb | df$tail>intervals_df$ub
  out_df = data.frame(p =df$p[out_idx],tail =df$tail[out_idx])
  
  # plot
  if (is.null(desc)){title = paste("Result for", n_distinct(df$boot) ,"samples")}
  else{title = desc}
  
  
  gg = ggplot(intervals_df) +
    geom_line(aes(x = p, y = mean, linetype = "Mean"), color = "black") +
    geom_line(aes(x = p, y = Q2, linetype = "Median"), color = "black") +
    # geom_ribbon(aes(x = p, ymin = Q1, ymax = Q3, fill = "1st and 3rd Quartiles"), alpha = 0.3) +
    geom_line(data = df,aes(x = p, y = tail,group = boot),col = 'grey', alpha = 0.3) +
    geom_line(data=original,aes(x=p,y=tail,color="Original"))+
    geom_point(data=original,aes(x=p,y=tail,color="Original"))+
    scale_x_log10() +
    labs(
      y = "tail index",
      x = "sample fraction for quantile estimation",
      title = title,
    ) +
    scale_linetype_manual(
      name = "Legend",
      values = c("Mean" = "dashed", "Median" = "solid"),
      labels = c("Mean", "Median")
    ) +
    scale_color_manual(name = "Legend",values=c("Original" = "red"), labels=c("Original"="Original dataset"))
  if(is.numeric(fix_y)){
    gg = gg + ylim(fix_y)
  }
  return(gg)
}

tailindexplot_5 <- function(boot,pconf=0.95,desc = NULL,fix_y=NULL){
  z <- abs(qnorm((1-pconf)/2))
  
  
  # print(true)
  # Calculate quartiles for each value of p
  intervals_df <- boot$df %>%
    group_by(p) %>%
    summarise(
      std = sd(tail)
    ) %>%
    ungroup()
  std = -cummax(-intervals_df$std)
  
  intervals_df[c('lb','ub')] <- cbind(boot$original$original$tail-z*std,boot$original$original$tail+z*std)
  print(intervals_df$p)
  
  # plot
  if (is.null(desc)){title = paste("Result for", n_distinct(boot$df$boot) ,"samples")}
  else{title = desc}

  
  gg = tailindexplot_2(list(boot$original))$gg +
    geom_ribbon(data = intervals_df,aes(ymin = 0, ymax = 1, color = "bootstrap"), alpha = 0.2)
  return(gg)
}


# stability plot of all columns from the dataset. (deprecated, rather use stabilityplot_col)
# timesep : time separation fr the peak selection (usually 12/timestep of the dataset)
# method : distribution "GP" or "GW"
tailindexplot_all_heights <- function(data,timesep,fraq_min=2.5,method='GP',fix_y = FALSE,plot_bm = FALSE){
    alpha = 0.05
    labels = names(data)#c('measurements','KNW','RACMO','RACMO','RACMO','RACMO')
    
    names(labels) = names(data)
    out = list()
    bm_fit = list()
    X = list()
    lower_bounds = c()
    upper_bounds =  c()
    
    for (i in names(data)){
      print(i)
      if (!is.null(data[[i]]) & !(i %in% c('DateTime','Year',grep('D[0-9]+|^[Puv]',names(data),value=TRUE))) || i=='PF010'){
        X[[i]] <-  PoTselect_2(na.approx(data[[i]]),0.1,timesep)
        l0 <- round(10^(seq(-fraq_min, -0, 0.05))*length(X[[i]]$pot))-1
        if (method == 'GP'){
          out[[i]] <- FitGP_MLE2(X[[i]]$pot, 1, N= 0, r11= 1, l0= l0, sigma= Inf, metadata= data.frame(label=i))
        }
        else if (method == 'GW'){
          out[[i]] <- FitGW_iHilli(X[[i]]$pot, 1, N= 0, r11= 1, l0= l0, sigma= Inf, metadata= data.frame(label=i))
        }
        else{
          cat("Error : no method called", blue(method),'\n')
          flush.console()
          return(NULL)
        }
      }
    }
    tip = tailindexplot_2(out,fix_y = TRUE,parameter = 'all')
    
    
    for (par in c('tail','scale','loc')){
      # plot
      gg = ggplot(tip[[par]]$df, aes(x = p, y = val, color = group)) +
        geom_point() +
        geom_line() +
        geom_ribbon(aes(ymin = lb, ymax = ub), alpha=0.02) +
        labs(x ="sample fraction for quantile estimation", y =par, color = "group") +
        scale_x_log10()+
        facet_wrap(~group,ncol=2, scales = "free_y")
      
      if (plot_bm){
        bm_df <- data.frame(group=unique(tip[[par]]$df$group))
        for (i in unique(tip[[par]]$df$group)){
          bm_fit[[i]] <- bootstrap(data,Nb=100,column = i,method='BM')
          bm_df[bm_df$group==i,'bm'] <- bm_fit[[i]]$original[par]
          bm_df[bm_df$group==i,'bm_lb'] <- bm_fit[[i]]$original[par] - qnorm(0.975)*sd(bm_fit[[i]]$df[[par]])
          bm_df[bm_df$group==i,'bm_ub'] <- bm_fit[[i]]$original[par] + qnorm(0.975)*sd(bm_fit[[i]]$df[[par]])
          bm_df[bm_df$group==i,'threshold'] <- sum(X[[i]]$pot > bm_fit[[i]]$original['loc'],na.rm=T)/length(data[[i]])
        }
        gg= gg +
        geom_vline(data=bm_df,aes(color=group,xintercept = threshold))+
        geom_hline(data=bm_df,aes(yintercept = bm,color=group))+
        geom_hline(data=bm_df,aes(yintercept = bm_lb,color=group))+
        geom_hline(data=bm_df,aes(yintercept = bm_ub,color=group))
      }
      
      if(is.numeric(fix_y)){
        gg = gg + ylim(fix_y)
      }
      else if (fix_y){
        scale_view_offset = (max(tip[[par]]$df$ub) - min(tip[[par]]$df$lb))/100
        gg = gg + ylim(min(tip[[par]]$df$lb)-scale_view_offset,max(tip[[par]]$df$ub)+scale_view_offset)
      }
      print(gg)
    }
}


# stability of data for all specified columns
#   data : whole dataset
#   timesep : time separation fr the peak selection (usually 12h/timestep of the dataset)
#   method : distribution "GP" or "GW"
#   columns : c(), columns you want to use for estimation
#   parameter : "tail" index, "scale", "location", or "all"
#   fraq_min : minimum sample fraction used for the estimations
#   th : if specified, the only output is the estimation for this threshold (as a sample fraction)  for all columns in `columns`
#   peak_frac: Percentage of the storm threshold used in PoT selection (default: 0.3)
#   winter: Boolean flag to consider only winter months (default: FALSE)
#   ncol : number of plot columns for the figure

stabilityplot_col <- function(data,timesep,fraq_min=2.5,columns = NULL,method='GP',fix_y = T,parameter='all',peak_frac=0.3,th=NULL,winter=F,ncol=3){
  alpha = 0.05
  
  out = list()
  bm_fit = list()
  X = list()
  lower_bounds = c()
  upper_bounds =  c()

  
  if (length(peak_frac)==1){peak_frac = rep(peak_frac,length(columns))}

  for (i in 1:length(columns)){
    print(paste(columns[i],peak_frac[i]))
    
    ws <- na.approx(data[[columns[i]]])
    if (winter){
      ws <- na.approx(data[[columns[i]]][month(data$DateTime) %in% c(10:12,1:3)])
    }
    print(length(ws))
    
    if (!is.null(ws)){
      X[[i]] <-  PoTselect_2(ws,peak_frac[i],timesep)
      if (is.null(th)) {l0 <- round(10^(seq(-fraq_min, -0, 0.05))*length(X[[i]]$pot))-1}
      else {l0 <- round(th*length(X[[i]]$pot))-1}
      if (method == 'GP'){
        out[[i]] <- FitGP_MLE2(X[[i]]$pot, 1, N= 0, r11= 1, l0= l0, sigma= Inf, metadata= data.frame(label=paste(columns[i],peak_frac[i])))
      }
      else if (method == 'GW'){
        out[[i]] <- FitGW_iHilli(X[[i]]$pot, 1, N= 0, r11= 1, l0= l0, sigma= Inf, metadata= data.frame(label=paste(columns[i],peak_frac[i])))
      }
      else{
        cat("Error : no method called", blue(method),'\n')
        flush.console()
        return(NULL)
      }
    }
  }

  tailindexplot_2(out,fix_y = fix_y,parameter = parameter,ncol = ncol)

}



# plot the correlation of the selected peaks for different sample fraction for the convergence threshold.
corrplot_col <- function(data,timesep,fraq_min=2.5,columns = NULL,peak_frac=0.3){
  alpha = 0.05
  
  out = list()
  bm_fit = list()
  X = list()
  lower_bounds = c()
  upper_bounds =  c()
  if (length(peak_frac)==1){peak_frac = rep(peak_frac,length(columns))}
  
  for (i in 1:length(columns)){
    print(paste(columns[i],peak_frac[i]))
    if (!is.null(data[[columns[i]]])){
      X[[i]] <-  PoTselect_2(na.approx(data[[columns[i]]]),peak_frac[i],timesep)
      l0 <- round(10^(seq(-fraq_min, -0, 0.05))*length(X[[i]]$pot))-1
      pval = c()
      for (l in l0){
        Xord = X[[i]]$pot[which(X[[i]]$pot >= X[[i]]$pot[l])]
        pval <- c(pval,Box.test(Xord , type="Ljung-Box")$p.value)
      }
      out[[i]] = list(l0=l0, pval= pval)
      }
  }

  par(mfrow = c(2, 2))  # Adjust the number of rows and columns as needed
  
  for (i in seq_along(out)) {
    plot(out[[i]]$l0, out[[i]]$pval, main = paste(columns[i],peak_frac[i]))
  }
  }




# plot the stability of the estimation for all parameters and 2 different datasets for the GP distribution
# 
stability_est_GP <- function(data_model, data_measure, col_model='w10m',col_measure='F010', timestep_model=3,timestep_measure=1/6,fix_y = FALSE,fraq_min=2, parameter='tail',overlap=F,peak_frac=0.3,winter = F,selected = 0.1,label=NULL){
  
  if(is.null(label)){
    label = c("model","measure","measure + fixed tail")
  }
  
  ws_model <- na.approx(data_model[[col_model]])
  ws_measure <- na.approx(data_measure[[col_measure]])
  if (winter){
    ws_model <- na.approx(data_model[[col_model]][month(data_model$DateTime) %in% c(10:12,1:3)])
    ws_measure <- na.approx(data_measure[[col_measure]][month(data_measure$DateTime) %in% c(10:12,1:3)])
  }
  
  pot_model <- PoTselect_2(ws_model,peak_frac,12/timestep_model)
  l0_model <- round(10^(seq(-fraq_min,0,0.05))*length(pot_model$pot))-1
  fit_model <- FitGP_MLE2(pot_model$pot, 1, N= 0, r11= 1, l0= l0_model,metadata = list(label=label[1]))
  
  
  selected = tailindexplot_2(list(fit_model),selected_threshold = selected,parameter = 'tail')
  
  # Function to constrain the tail index parameter in the estimation
  fixedpar = list(gamma0 = selected$selected['tail'],gamma0Std = selected$selected['sd'])
  
  
  pot_measure <- PoTselect_2(ws_measure,peak_frac,12/timestep_measure)
  l0_measure <- round(10^(seq(-fraq_min,0,0.05))*length(pot_measure$pot))-1
  fit_measure <- FitGP_MLE2(pot_measure$pot, 1, N= 0, r11= 1, l0= l0_measure,metadata = list(label=label[2]))
  
  if (!is.null(label) && length(label)==2){plot_list = list(fit_model,fit_measure)}
  else {
    fit_measure_model <- FitGP_MLE2(pot_measure$pot, 1, N= 0,fixedpar = fixedpar, r11= 1, l0= l0_measure,metadata = list(label=label[3]))
    plot_list = list(fit_model,fit_measure,fit_measure_model)
    }
  tailindexplot_2(plot_list,ncol=2,fix_y = fix_y,parameter = parameter,overlap = overlap)
}





# plot the stability of the estimation for all parameters and 2 different datasets for the GW distribution

stability_est_GW <- function(data_model, data_measure, col_model='w10m',col_measure='F010', timestep_model=3,timestep_measure=1/6,fix_y = FALSE, fraq_min=2, parameter='tail',overlap=F,peak_frac=0.3,winter=F,selected = 0.1,label=NULL){
  if(is.null(label)){
    label = c("model","measure","measure + fixed tail")
  }
  
  ws_model <- na.approx(data_model[[col_model]])
  ws_measure <- na.approx(data_measure[[col_measure]])
  if (winter){
    ws_model <- na.approx(data_model[[col_model]][month(data_model$DateTime) %in% c(10:12,1:3)])
    ws_measure <- na.approx(data_measure[[col_measure]][month(data_measure$DateTime) %in% c(10:12,1:3)])
  }
  
  
  pot_model <- PoTselect_2(ws_model,peak_frac,12/timestep_model)
  l0_model <- round(10^(seq(-fraq_min,0,0.05))*length(pot_model$pot))-1
  fit_model <- FitGW_iHilli(pot_model$pot, 1, N= 0, r11= 1, l0= l0_model,metadata = list(label=label[1]))
  
  
  selected = tailindexplot_2(list(fit_model),selected_threshold = selected,parameter = 'tail')

  # Function to constrain the tail index parameter in the estimation
  fixedpar = list(theta0 = selected$selected['tail'],theta0Std = selected$selected['sd'])
  
  
  pot_measure <- PoTselect_2(ws_measure,peak_frac,12/timestep_measure)
  l0_measure <- round(10^(seq(-fraq_min,0,0.05))*length(pot_measure$pot))-1
  fit_measure <- FitGW_iHilli(pot_measure$pot, 1, N= 0, r11= 1, l0= l0_measure,metadata = list(label=label[2]))
  
  if (!is.null(label) && length(label)==2){plot_list = list(fit_model,fit_measure)}
  else {
    fit_measure_model <- FitGW_iHilli(pot_measure$pot, 1, N= 0,fixedpar = fixedpar, r11= 1, l0= l0_measure,metadata = list(label=label[3]))
    plot_list = list(fit_model,fit_measure,fit_measure_model)
    }
  tailindexplot_2(plot_list,ncol=2,fix_y = fix_y,parameter = parameter,overlap = overlap)
}