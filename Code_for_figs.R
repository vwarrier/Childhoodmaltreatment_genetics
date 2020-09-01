###Install packages
library("CMplot")
library("ggplot2")
library(ggcorrplot)


## Fig 2

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # colour blind colour palette

a = ggplot(data1, aes(x = Phenotype, y = h2SNP, fill = Phenotype, colour = Phenotype)) + 
  geom_point(size = 5) + geom_errorbar(width=.3, aes(ymin = h2SNP-1.96*SE, ymax = h2SNP+1.96*SE)) +
  theme_classic() + ylab("SNP heritability") + geom_hline(yintercept = 0) + xlab("Phenotype") +   scale_colour_manual(values=cbbPalette) +
  theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())

b = ggcorrplot(data1, outline.col = "white", lab = TRUE, colors = c( "#56B4E9", "white", "#E69F00"))

c = ggplot(data1, aes(x = Phenotype, y = h2SNP, fill = Phenotype, colour = Phenotype)) + 
  geom_point(size = 5) + geom_errorbar(width=.3, aes(ymin = h2SNP-1.96*SE, ymax = h2SNP+1.96*SE)) +
  theme_classic() + ylab("SNP heritability") + geom_hline(yintercept = 0) + xlab("Phenotype") +   scale_colour_manual(values=cbbPalette) +
  theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())

d = ggcorrplot(data1, outline.col = "white", lab = TRUE, colors = c( "#56B4E9", "white", "#E69F00"), hc.order = TRUE)


multiplot(a, c, b, d, cols=2) # See below for multiplot key


### Figure 3


CMplot(prosp,type="p",plot.type="q",r=0.4,col=c("grey30","grey60"),chr.labels=paste("Chr",c(1:22),sep=""),
       threshold=c(5e-8),cir.chr.h=1.5,amplify=FALSE,threshold.lty=c(1,2),threshold.col="red",signal.line=1,signal.col="red",
       bin.size=1e6,outward=TRUE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=8,height=8)

retro_prosp = retro_prosp[,c("SNP", "CHR", "BP", "P")]


#### Supplementary Figure - PGS 


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # colour blind friendly colour palette


pd <- position_dodge(width = 0.5)

ggplot(data, aes(x = Time, y = OR, fill = Model, colour = Model)) + 
  geom_point(size = 5, position = pd) + geom_errorbar(width=.3, aes(ymin = LCI, ymax = UCI), position = pd) +
  theme_classic() + ylab("ORs") + geom_hline(yintercept = 1) + xlab("Age range") +   scale_colour_manual(values=cbbPalette) 


#### Multiplot
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}