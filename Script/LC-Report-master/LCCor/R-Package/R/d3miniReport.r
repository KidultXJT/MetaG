d3hm <- function(
  df,
  brwCol = "RdYlBu",  # select From RColorBrewer::Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd 
  # Purples RdPu Reds YlGn YlGnBu YlOrBr YlOrRd BrBG PiYG PRGn PuOr
  # RdBu RdGy RdYlBu RdYlGn Spectral
  scale = "none",
  nolabCol = T,
  nolabRow = T,
  dendrogram = "none"
  ){
  # Package
  require(d3heatmap)
  
  if (nolabCol){labCol = rep("",ncol(df))}
  if (nolabCol){labRow = rep("",nrow(df))}
  return(d3heatmap(df,
            scale = scale,
            colors = colorRampPalette(rev(brewer.pal(n = 7, name = brwCol)))(100),
            xaxis_height= 120,
            yaxis_width = 120,
            dendrogram = dendrogram
            ))
  
}