ggHeatmap <- function(
  matrix,
  brwCol = "RdYlBu",  # select From RColorBrewer::Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd 
  # Purples RdPu Reds YlGn YlGnBu YlOrBr YlOrRd BrBG PiYG PRGn PuOr
  # RdBu RdGy RdYlBu RdYlGn Spectral
  color = "white",
  xtitle = "",
  ytitle = "",
  filltitle = "",
  nbin = 100,
  coord_flip = F,
  numDisplay = T,
  titlesize = 10,
  textsize  = 6,
  numsize = 2
){
  # Package
  require(ggplot2)
  
  df <- data.frame(Name=row.names(matrix),stack(matrix))
  df$Name <- factor(df$Name, levels=rev(unique(df$Name)))
  df$ind <- factor(df$ind, levels=unique(df$ind))
  
  p = ggplot(df, aes(x=factor(Name),y=ind,fill=values,label=round(values,2))) +
    geom_tile(color=color,size=.3) +
    labs(x=xtitle,y=ytitle,fill=filltitle) + 
    scale_fill_distiller(palette = brwCol,
                         guide = guide_colorbar(nbin=nbin,draw.ulim = FALSE,draw.llim = FALSE)) + reportStyle()
  
  if(numDisplay){p <- p + ggBasicText(Col="white",size = numsize)}
  # Text Size
  theme = theme(axis.title  = element_text(size=titlesize,colour="grey10",hjust = 0.5),
                axis.text.x = element_text(size=textsize,colour="grey40",hjust = 1,angle=90),
                axis.text.y = element_text(size=textsize,colour="grey40",hjust = 1),
                axis.line   = element_blank(),
                axis.ticks  = element_blank(),
                # BackGround
                panel.background = element_rect(fill = "white", color="white"))
  p = p + theme
  if(coord_flip){p <- p + coord_flip()}
  
  return(p)
}
