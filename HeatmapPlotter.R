library(ggdendro)
library(plyr)
library(ggplot2)
library(grid)
library(gtable)
library(colorspace)
library(RColorBrewer)
library(reshape2)
library(gridExtra)
library(cowplot)
suppressMessages(library("ggplot2"))
options("width"=180) # useful for debugging

# Return the command line argument associated with a given flag (i.e., -o foo),
# or the default value if argument not specified.
# Note that this will break if an argument is not supplied after the flag.
get_val_arg = function(args, flag, default) {
    ix = pmatch(flag, args)
    if (!is.na(ix)){ val = args[ix+1] } else { val = default }
    return(val)
}

# Return boolean specifying whether given flag appears in command line (i.e., -o),
get_bool_arg = function(args, flag) {
    ix = pmatch(flag, args)
    if (!is.na(ix)){ val = TRUE } else { val = FALSE }
    return(val)
}

# Usage: 
#   args = parse_args()
#   print(args$disease.filter)
parse_args = function() {
    args = commandArgs(trailingOnly = TRUE)

    # optional arguments
    verbose = get_bool_arg(args, "-v")
    height = as.numeric(get_val_arg(args, "-h", 8)) 
    width = as.numeric(get_val_arg(args, "-w", 6))  
    skip.x.label = get_bool_arg(args, "-x")
    skip.y.label = get_bool_arg(args, "-y")
    x.title = get_val_arg(args, "-X", "Sample")
    y.title = get_val_arg(args, "-Y", "Gene")
    gene.grouping = get_val_arg(args, "-g", NA)     
    sampleset.grouping = get_val_arg(args, "-p", NA)
    plot.title = get_val_arg(args, "-t", "Plot")    
    hide.legend = get_bool_arg(args, "-l")
    cluster.analysis = get_bool_arg(args, "-c")
    gradiant.scale = get_bool_arg(args, "-s")
    
    # mandatory positional arguments.  These are popped off the back of the array, last one listed first.
    plot.fn = args[length(args)];               args = args[-length(args)]
    data.fn = args[length(args)]; args = args[-length(args)]
    
    
    val = list( 'plot.fn'=plot.fn, 'verbose'=verbose, 'data.fn'=data.fn,
                'height'=height, 'width'=width, 'skip.x.label'=skip.x.label, 'skip.y.label'=skip.y.label,
                'gene.grouping'=gene.grouping, 'sampleset.grouping'=sampleset.grouping,
                'plot.title'=plot.title, 'hide.legend'=hide.legend, 'cluster.analysis'=cluster.analysis,
                'gradiant.scale'=gradiant.scale)
    
    if (val$verbose) { print(val) }
    return (val)
}

options("width"=180) # useful for debugging
#options("width"=270) 
args = parse_args()

# Read input file(s)
data<-read.table(args$data.fn, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, comment.char = "#") #Reads TSV
data.original = data
if(!is.na(args$gene.grouping)){
  gene_categories <- read.table(as.character(args$gene.grouping), header = FALSE, sep = "\t", comment.char = "#")
}
if(!is.na(args$sampleset.grouping)){
  sample_categories <- read.table(as.character(args$sampleset.grouping), header=FALSE, sep = "\t", comment.char = "#")
}

#PlotMaster H.A.K. [Hugs And Kisses] Felix Hu, Justin Chen ---------------------------------------------

#Functions

#Returns a palette of size (n) 
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

#Split string by ';' and returns either first or second substring
labSplit<-function(str, front = TRUE){
  t<-unlist(strsplit(str,";"))
  if(front) {
    return(t[1])
  } else {
    return(t[2])
  }
}

#Generate plots to act as side bars
bars <- function(dat, RColBrewPalette = "Set1", isTop = TRUE){

  getPalette = colorRampPalette(brewer.pal(brewer.pal.info[RColBrewPalette,]$maxcolors, RColBrewPalette))
  
  if(isTop){
    
    val <- dat[apply(dat, 1, function(r) any(r %in% as.character(dat$Gene[1]))),]
    f <- ggplot(val, aes(val, x = Sample, y = Gene)) + geom_tile(aes(fill = as.character(Sample_Category))) + scale_fill_manual(values = c(getPalette(length(unique(val$Sample_Category)))), name = "Gene Key")  + ggtitle(args$plot.title)
    
  } else {
    
    val <- dat[apply(dat, 1, function(r) any(r %in% as.character(dat$Sample[1]))),]
    f <- ggplot(val, aes(val, x = Sample, y = Gene)) + geom_tile(aes(fill = as.character(Gene_Category))) + scale_fill_manual(values = c(getPalette(length(unique(val$Gene_Category)))), name = "Sample Key")
    
  }
  
  f <- f + theme_bw()
  f <- f + theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), panel.border = element_blank(), legend.direction = "horizontal", plot.margin = unit(c(0, 0.03, 0, 0), "in"))
  
  return(f)
}

#Returns indices to sort data by method specified in the argument.
indices<-function(dat, method = "input", xy = 1){    #sorts `data`,  xy = 1 sorts x axis, xy = 2 sorts y axis
  #Makes a copy of dat 
  dataCl = dat                                     
  
  if(method == "alphabet"){                      
    if(xy == 1){
      val <- as.data.frame(as.character(colnames(dataCl)))
      names(val) = c("Sample")
      val$ind <- c(1:length(val$Sample))
      val <- arrange(val, Sample)
      return (val$ind)
    } else {
      val <- as.data.frame(as.character(rownames(dataCl)))
      names(val) = c("Gene")
      val$ind <- c(1:length(val$Gene))
      val <- arrange(val, Gene)
      return (val$ind)
    }
  } else if(method == "input") {
    if(xy == 1) { 
      return(c(1:ncol(dataCl))) 
    } else { 
      return(c(1:nrow(dataCl))) 
    }
  } else if(method == "clust") {
    if(xy == 1) {
      return (hclust(dist(t(dataCl)))$order)
    } else {
      return (hclust(dist(dataCl))$order)
    }
  } else if(method == "alphabetgrouping"){
    if(xy == 1 && !is.na(args$sampleset.grouping)){
      val <- as.data.frame(as.character(colnames(dataCl)))
      names(val) = c("Sample")
      val$ind <- c(1:nrow(val))
      names(sample_categories) = c("Sample", "Sample_Category")
      val <- merge(x = val, y = sample_categories, all.x = TRUE)
      val$interac <- interaction(val$Sample_Category, val$Sample, lex.order = FALSE, sep = ";")
      val <- arrange(val, as.character(interac))
      return(val$ind)
    } else if(!is.na(args$gene.grouping)){
      val <- as.data.frame(as.character(rownames(dataCl)))
      names(val) = c("Gene")
      val$ind <- c(1:nrow(val))
      names(gene_categories) = c("Gene", "Gene_Category")
      val <- merge(x = val, y = gene_categories, all.x = TRUE)
      val$interac <- interaction(val$Gene_Category, val$Gene, lex.order = FALSE, sep = ";")
      val <- arrange(val, as.character(interac))
      return(val$ind)
    } else if(xy == 1){
      return(c(1:ncol(data)))
    } else {
      return(c(1:nrow(data)))
    }
  }
}

# Function for panel alignment from Baptiste's gridExtra library.
# Found at https://github.com/baptiste/gridextra/wiki/arranging-ggplot
gtable_frame <- function(g, width=unit(1,"null"), height=unit(1,"null")){
  panels <- g[["layout"]][grepl("panel", g[["layout"]][["name"]]), ]
  ll <- unique(panels$l)
  tt <- unique(panels$t)
  
  fixed_ar <- g$respect
  if(fixed_ar) { # Align with aspect ratio constraints
    ar <- as.numeric(g$heights[tt[1]]) / as.numeric(g$widths[ll[1]])
    print(ar)
    height <- width * ar
    g$respect <- FALSE
  }
  
  core <- g[seq(min(tt), max(tt)), seq(min(ll), max(ll))]
  top <- g[seq(1, min(tt)-1), ]
  bottom <- g[seq(max(tt)+1, nrow(g)), ]
  left <- g[, seq(1, min(ll)-1)]
  right <- g[, seq(max(ll)+1, ncol(g))]
  
  fg <- nullGrob()
  lg <-  if(length(left))  g[seq(min(tt), max(tt)), seq(1, min(ll)-1)] else fg
  rg <- if(length(right)) g[seq(min(tt), max(tt)), seq(max(ll)+1,ncol(g))] else fg
  grobs = list(fg, g[seq(1, min(tt)-1), seq(min(ll), max(ll))], fg, 
               lg, g[seq(min(tt), max(tt)), seq(min(ll), max(ll))], rg, 
               fg, g[seq(max(tt)+1, nrow(g)), seq(min(ll), max(ll))], fg)
  widths <- unit.c(sum(left$widths), width, sum(right$widths))
  heights <- unit.c(sum(top$heights), height, sum(bottom$heights))
  all <- gtable_matrix("all", grobs = matrix(grobs, ncol=3, nrow=3, byrow = TRUE), 
                       widths = widths, heights = heights)
  all[["layout"]][5,"name"] <- "panel" # make sure knows where the panel is
  if(fixed_ar)  all$respect <- TRUE
  all
}

create.dendrogram <- function(dat, xy = 1) { #xy = 1 operates on x axis, xy = 2 operates on y axis
  if(xy == 1) {
   dd <- dendro_data(as.dendrogram(hclust(dist(t(dat)))))
   dendrogram.plot <- ggplot(segment(dd)) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend), color="gray50") + scale_y_reverse() + theme_dendro() + theme(plot.margin=unit(c(-0.05,-0.0475,-0.2,-0.0505), "npc"))
  } else {
    dd <- dendro_data(as.dendrogram(hclust(dist(dat))))
    dendrogram.plot <- ggplot(segment(dd)) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend), color="gray50") +  coord_flip() + scale_y_reverse() + theme_dendro() + theme(plot.margin=unit(c(0,-0.1,0,0), "npc"))
  }
  return(dendrogram.plot)
}
  
met = "clust"

#Reformat 'data' samples
data <- data[,indices(data, method = met)]

#Reformat 'data' genes
data <- data [indices(data, xy = 2, method = met),]

#Read plot data
plotData<-as.data.frame(melt(as.matrix(data)))
names(plotData) = c('Gene', 'Sample', 'Variant')
#Keep track of original indices of plotData
plotData$ind <- c(1:nrow(plotData))

#Gene grouping
if(!is.na(args$gene.grouping)){
  names(gene_categories) = c("Gene", "Gene_Category")
  
    #Error checking
    if(length(gene_categories$Gene) < length(unique(plotData$Gene))){
      stop("Your plotting data and grouping manual do not match: \"", c(setdiff(as.character(gene_categories$Gene), as.character(unique(plotData$Gene))), setdiff(as.character(unique(plotData$Gene)), as.character(gene_categories$Gene))), "\" fail to overlap")
    }
  
  plotData<-merge(x = plotData, y = gene_categories, all.x = TRUE)
}

#Sampleset grouping
if(!is.na(args$sampleset.grouping)){
  names(sample_categories) = c("Sample", "Sample_Category")
  
    #Error checking
    if(length(sample_categories$Sample) < length(unique(plotData$Sample))){
      stop("Your plotting data and sampleset grouping manual do not match: \"", c(setdiff(as.character(sample_categories$Sample), as.character(unique(plotData$Sample))), setdiff(as.character(unique(plotData$Sample)), as.character(sample_categories$Sample))), "\" fail to overlap")
    }  
  
  plotData<-merge(x = plotData, y = sample_categories, all.x = TRUE)
}

#Sort plotData according to original indices
plotData <- plotData[order(plotData$ind),]

#Prepare y-axis labels if gene grouping is true
if(!is.na(args$gene.grouping)){
  labY = unique(sort(interaction(plotData$Gene_Category, plotData$Gene, lex.order = FALSE, sep = ";")))
  lY2<-unlist(lapply(labY, function(x) labSplit(as.character(x))))
  lY<-unlist(lapply(labY, function(x) labSplit(as.character(x), FALSE)))
}

#Prepare x-axis labels if sampleset grouping is true
if(!is.na(args$sampleset.grouping)){
  labX = unique(sort(interaction(plotData$Sample_Category, plotData$Sample, lex.order = FALSE, sep = ";")))
  lX2<-unlist(lapply(labX, function(x) labSplit(as.character(x))))
  lX<-unlist(lapply(labX, function(x) labSplit(as.character(x),front = FALSE)))
}

#Heatmap Plot
#Tell ggPlot that plotData$Sample is a factor
plotData$Sample <- as.character(plotData$Sample)
plotData$Sample <- factor(plotData$Sample, levels = unique(plotData$Sample))

#Tell ggPlot that plotData$Gene is a factor
plotData$Gene <- as.character(plotData$Gene)
plotData$Gene <- factor(plotData$Gene, levels = unique(plotData$Gene))

p <- ggplot(plotData, aes(x=Sample, y=Gene))

if(args$gradiant.scale){
  p = p + geom_tile(aes(fill = Variant)) + scale_fill_gradient(low = "white", high = "black", name = "Key") + theme_bw()
} else {
  p = p + geom_tile(aes(fill = as.character((Variant)))) + scale_fill_manual(values = c("white", getPalette(length(unique(plotData$Variant)))), name = "Key") + theme_bw()
}

#Format x axis labels
if(args$skip.x.label == TRUE){
  textx = element_blank()
} else {
  textx = element_text(angle=-90, hjust=0, vjust=1, color = "black")
}

#Format y axis labels
if(args$skip.y.label == TRUE){
  texty = element_blank()
} else {
  texty = element_text(color = "black")
}

p = p + theme(axis.title.x = element_blank(), axis.text.x = textx, axis.title.y = element_blank(), axis.text.y = texty, panel.grid = element_blank(),
        axis.ticks = element_blank(),  panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        plot.margin = unit(c(0.01, 0.01, 0, 0),"npc"), legend.direction = "horizontal", legend.position = "right")

#Legend
if(args$hide.legend == TRUE){
  p = p + guides(fill = FALSE)
}
  
#Top bar sampleset grouping
if(!is.na(args$sampleset.grouping)){
  
  #Get top bar coordinates
  indX<-c(which(diff(as.numeric(as.factor(lX2)))!=0))
  indX = c(indX, length(lX2))

}

#Side bar gene grouping
if(!is.na(args$gene.grouping)){
  
  #Get side bar coordinates
  indY<-c(which(diff(as.numeric(as.factor(lY2)))!=0))
  indY = c(indY, length(lY2))
  
}

#Build heatmap plot
g <- ggplot_gtable(ggplot_build(p + guides(fill = FALSE)))

#g <- gtable_add_rows(g, unit(5,"cm"))
#g <- gtable_add_grob(g, ggplotGrob(create.dendrogram(data.original, xy = 1)), t = nrow(g), l=4, b=nrow(g), r=4)
#g <- gtable_add_cols(g, unit(5, "cm"), pos = 0)

g2 <- ggplot_gtable(ggplot_build(p))
g$layout$clip[g$layout$name == "panel"] <- "off"
legend1 <- gtable_filter(g2, "guide-box")

print(plotData)

#Build topbar plot
l1 <- bars(dat = plotData, RColBrewPalette = "Set2", isTop = T)
s1 = l1 + guides(fill = FALSE)
d1 <- ggplot_gtable(ggplot_build(s1))
k1 <- ggplot_gtable(ggplot_build(l1))
d1$layout$clip[d1$layout$name == "panel"] <- "off"
legend2 <- gtable_filter(k1, "guide-box")
d1$widths = g$widths

#Build sidebar plot
l2 <- bars(dat = plotData, isTop = F) 
s2 = l2 + guides(fill = FALSE)
k2 <- ggplot_gtable(ggplot_build(l2))
d2 <- ggplot_gtable(ggplot_build(s2))
d2$layout$clip[d2$layout$name == "panel"] <- "off"
legend3 <- gtable_filter(k2, "guide-box")
d2$heights = g$heights

#Combine plots
blank <- ggplot_gtable(ggplot_build(ggplot() + geom_blank())) # Blank ggplot
blank$heights = d1$heights
blank$widths = d2$widths

#Combine 3 legends
legendGrob <- gtable:::cbind_gtable(legend1, legend2, "last")
legendGrob2 <- gtable:::cbind_gtable(legendGrob, legend3, "last")

fg1 <- gtable_frame(d1, width = unit(31,"null"), height = unit(1,"null"))
fg2 <- gtable_frame(d2, width = unit(1, "null"), height = unit(31, "null"))
fg3 <- gtable_frame(g, width = unit(31, "null"), height = unit(31, "null"))

fgblank <- gtable_frame(blank, width = unit(10, "null"), height = unit(1, "null"))

fg13 <- gtable_frame(rbind(fg1, fg3), width=unit(31,"null"), height=unit(32,"null"))

fg2blank <- gtable_frame(rbind(fgblank, fg2), width = unit(1,"null"), height = unit(32,"null"))

combined.group <- gtable_frame(cbind(fg13, fg2blank), width = unit(32, "null"), height = unit(32, "null"))



#------------------------------------------------------------------------------------------------------------
# Save
cat(sprintf("Saved to %s\n", args$plot.fn))
pdf(args$plot.fn, height=args$height, width=args$width)

#gtable_show_layout(g)
#grid::grid.draw(k1)

#grid::grid.draw(create.dendrogram(data, xy = 1))
#grid::grid.draw(combined.group)
#grid.arrange(combined.group, bottom.dendrogram.gt, side.dendrogram.gt, widths = c(2, 31, 1), heights = c(1, 31, 2), layout_matrix = rbind(c(NA,1,1),
#                                                                                                                                        c(3,1,1),
#                                                                                                                                        c(NA,2,NA)))

#Closer answers
#main.plot.gb <- arrangeGrob(d1, blank, g, d2, widths = c(45, 1), heights = c(1, 20), ncol = 2, nrow = 2) #layout_matrix = rbind(c(1, NA), c(1, 2)))

#main.plot.gb <- arrangeGrob(d1, g, d2, legendGrob2, widths = c(45, 1), heights = c(1, 20, 1), layout_matrix = rbind(c(1, NA),
#                                                                                                                    c(2, 3),
#                                                                                                                    c(4, 4)))


#plot_grid(d1, blank, g, d2, align = "hv", nrow = 2, ncol = 2, rel_widths = c(10, 1), rel_heights = c(1, 2))

#Closest answer
grid.arrange(combined.group, legendGrob2, heights = c(8, 1), nrow = 2, newpage = FALSE)

grid.newpage()
gtable_show_layout(combined.group)
dev.off()

