ggBasicPoint <- function(
  Col=NULL, # When Use One Color
  size=4,
  alpha=.7,
  df=NULL,  # When Use This, Try To make a Basic Graph; Else, To make a Point Graph LAYER !!! 
  x=1,
  y=2
){
  # Description:
  # A point Graph for Basic Graph or A Layer
  
  # Package:
  require(ggplot2)
  
  if(is.null(Col)){
    point = geom_point(size = size, alpha = alpha)
  }else{
    point = geom_point(size = size, alpha = alpha, color = Col)
  }
  
  if(!is.null(df)){
    point = ggplot(data = df,
                  aes(x=df[,x],
                      y=df[,y])) + point
  }
  
  return(point)
}

ggBasicText <- function(
  Col=NULL, # When Use One Color
  size=4,
  df=NULL,  # When Use This, Try To make a Basic Graph; Else, To make a Text Graph LAYER !!! 
  x=1,
  y=2,
  label=3
){
  # Description:
  # A TEXT Graph for Basic Graph or A Layer
  # Package:
  require(ggplot2)
  
  if(is.null(Col)){
    text = geom_text(size = size)
  }else{
    text = geom_text(size = size,
                     color = Col)
  }
  
  if(!is.null(df)){
    
    df$x = factor(df$x)
    text = ggplot(data = df,
                  aes(x=df[,x],
                      y=df[,y],
                      label=df[,label])) + 
    text
  }
  
  return(text)
}