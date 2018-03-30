## --------------------------- ##
##   Title: cor::Heatmap       ##
##  Author: Kidult             ##
##    Date: 2017/9/27          ##
## --------------------------- ##
reportStyle <- function(
  bgCol = "white"
  ){
  # Description:
  # Draw A Graph With Sagene Report Style. (By GGplot2::theme)
  # 
  # Package
  require(ggplot2)
  
  theme(# Panel
        panel.border  = element_blank(),
        panel.grid    = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = bgCol, colour = "grey20"),
        
        # Axis
        axis.title  = element_text(size=16,colour="grey10",hjust = 0.5),
        axis.text   = element_text(size=10,colour="grey20"),
        axis.line   = element_line("grey40"),
        
        # Legend
        legend.position     ="right",
        legend.justification="top",
        legend.title      = element_text(face = "bold"),
        legend.key.size   = unit(.4, "cm"),
        legend.key.height = unit(.4, "cm"),
        legend.key.width  = unit(.4, "cm"),
        legend.text      = element_text(angle=0,
                                        hjust=.95,
                                        vjust=.95,
                                        size=8,
                                        colour = "grey56"),
        
        # Strip
        strip.text       = element_text(size=14,colour="grey10",hjust = 0.5),
        strip.background = element_rect(colour="grey20", fill="grey90")
        )
}

colScale <- function(
  nColors, 
  UrColPallet = NULL # Example :: c("#001F3F","#FF5722","#51710A")
){
  # Description 
  # If You Use This you can Scale_color_* again !!!!
  # UrColPallet Should more than Your df[,Name]
  
  if(is.null(UrColPallet)){
    UrColPallet =  c("#001F3F","#FF5722","#51710A","#4E1184","#F87D09","#C5D200",
                     "#932B77","#0D63A5","#F9C535","#9BDF46","#FD367E","#228896",
                     "#F8DA5B","#83CC61","#E84A5F","#3D84A8","#F8C957","#B7E576",
                     "#FF847C","#46CDCF","#99CDA9","#FECEA8","#ABEDD8","#D1E9D2",
                     "#FFD6A4","#A7CDCC","#E5F4E7","#FDE9DF","#D6E6F2","#F1FDF3",
                     "#43496E","#544D7E","#65589C","#3D3551","#884EA2","#3EC280",
                     "#D15400","#EDEDEE","#E16A6B","#D9DEE1","#EA9432","#36D6B5",
                     "#68C2A1","#AAB6B6","#2474A8","#3597DA","#E57E22","#81CEE0",
                     "#D54542","#D24D58","#51B2D7","#F5AA36","#AEA7D1","#00B069",
                     "#EE4836","#E77E04","#F6C917","#1F3A92","#2ABA99","#169F85",
                     "#F5D66D","#BDC2C6","#F17835","#1E8AC2","#F44746","#A0DDCE",
                     "#66CB99","#59AAE2","#903D87","#D81E17","#8C44AC","#4ECCC3",
                     "#F22613","#DA0A5B","#9959B4","#95281B","#03C8A7","#F1F0EE",
                     "#BEBEBE","#F62459","#94A4A4","#C7F5C4","#D1D6D3","#F89306",
                     "#EA964E","#65C5BA","#5C96BE","#23A6F0")

  }else{
    UrColPallet = UrColPallet
  }
  
  UrColPallet = rep(UrColPallet,100)

  s = scale_color_manual(values = c(UrColPallet[1:nColors]))
  
  return(s)
}

fillScale <- function(
  nColors, 
  UrColPallet = NULL # Example :: c("#001F3F","#FF5722","#51710A")
){
  # Description 
  # If You Use This you can Scale_color_* again !!!!
  # UrColPallet Should more than Your df[,Name]
  
  if(is.null(UrColPallet)){
    UrColPallet =  c("#001F3F","#FF5722","#51710A","#4E1184","#F87D09","#C5D200",
                     "#932B77","#0D63A5","#F9C535","#9BDF46","#FD367E","#228896",
                     "#F8DA5B","#83CC61","#E84A5F","#3D84A8","#F8C957","#B7E576",
                     "#FF847C","#46CDCF","#99CDA9","#FECEA8","#ABEDD8","#D1E9D2",
                     "#FFD6A4","#A7CDCC","#E5F4E7","#FDE9DF","#D6E6F2","#F1FDF3",
                     "#43496E","#544D7E","#65589C","#3D3551","#884EA2","#3EC280",
                     "#D15400","#EDEDEE","#E16A6B","#D9DEE1","#EA9432","#36D6B5",
                     "#68C2A1","#AAB6B6","#2474A8","#3597DA","#E57E22","#81CEE0",
                     "#D54542","#D24D58","#51B2D7","#F5AA36","#AEA7D1","#00B069",
                     "#EE4836","#E77E04","#F6C917","#1F3A92","#2ABA99","#169F85",
                     "#F5D66D","#BDC2C6","#F17835","#1E8AC2","#F44746","#A0DDCE",
                     "#66CB99","#59AAE2","#903D87","#D81E17","#8C44AC","#4ECCC3",
                     "#F22613","#DA0A5B","#9959B4","#95281B","#03C8A7","#F1F0EE",
                     "#BEBEBE","#F62459","#94A4A4","#C7F5C4","#D1D6D3","#F89306",
                     "#EA964E","#65C5BA","#5C96BE","#23A6F0")
    
  }else{
    UrColPallet = UrColPallet
  }
  
  UrColPallet = rep(UrColPallet,100)
  
  s = scale_fill_manual(values = c(UrColPallet[1:nColors]))
  
  return(s)
}
