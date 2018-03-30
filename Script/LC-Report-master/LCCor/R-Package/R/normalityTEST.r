## ------------------------- ##
##   Title: Normality TEST   ##
##  Author: Kidult           ##
##    Date: 2017/11/13       ##
## ------------------------- ##
NormalTEST <- function(
  X, # A Table
  method = "shapiro" # shapiro (means W Test) ks (means D Test)
){
  # Description: 
  
  # 
  # Package:
  require(nortest)
  
  if(method == "shapiro"){
    wlst = c()
    for(x in colnames(X)){
      w = shapiro.test(X[,colnames(X)==x])$p.value
      wlst = c(wlst,w)
    }
    D = data.frame(Var = colnames(X),ShapiroTest = wlst)
  }else if(method == "ks"){
    dlst = c()
    for(x in colnames(X)){
      if(length(X[,colnames(X)==x]) > 4){
        d = lillie.test(X[,colnames(X)==x])$p.value # ks.test Should provide a Y vector
        dlst = c(dlst,d)
      }else{dlst = c(dlst,"NA")}
      }
    D = data.frame(Var = colnames(X),KSTest = dlst)
    }
  return(D)
}

QQ <- function(
  X # A Vector
){
  # Description: 
  # 
  # Make QQ Plot
  # 
  # Package:
  require(ggplot2)
  
  qqplot = ggplot() +
    geom_qq(aes(sample = X),
            color = 'black',
            size  = 2,
            shape = 1,
            alpha = .6) +
    geom_abline(intercept = mean(X),
                slope     = sd(X),
                color     = 'red',
                size      = .5,
                alpha     = .8) +
    ggtitle(label = 'Q-Q Plot') +
    xlab(label    = 'Sample Quantile') +
    ylab(label    = 'Theoretical Quantile') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          axis.title = element_text(size = 14),
          axis.text  = element_text(size = 12))
  
  return(qqplot)
}

PP <- function(
  X # A Vector
){
  # Description: 
  # 
  # Make PP Plot
  # 
  # Package:
  require(ggplot2)
  
  n = length(X)
  p = (1 : n) / n - 0.5 / n
  
  ppplot = ggplot() +
    geom_point(aes(y = p,
                   x = sort(pnorm(X, mean(X), sd(X)))),
               color = 'black',
               size = 2,
               shape = 1,
               alpha = .6) +
    geom_abline(slope = 1,
                color = 'red',
                size = .5,
                alpha = .8) +
    xlim(0,1) + ylim(0,1) +
    ggtitle(label = 'P-P Plot') +
    xlab(label = 'Sample Percentile') +
    ylab(label = 'Theoretical Percentile') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))
  
  return(ppplot)
}

