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
    hide.legend = get_bool_arg(args, "-l") #To implement
    gradiant.scale = get_bool_arg(args, "-s")
    sort.method = as.character(get_val_arg(args, "-m", "alphabetgrouping"))
    
    # mandatory positional arguments.  These are popped off the back of the array, last one listed first.
    plot.fn = args[length(args)];               args = args[-length(args)]
    data.fn = args[length(args)]; args = args[-length(args)]
    
    
    val = list( 'plot.fn'=plot.fn, 'verbose'=verbose, 'data.fn'=data.fn,
                'height'=height, 'width'=width, 'skip.x.label'=skip.x.label, 'skip.y.label'=skip.y.label,
                'gene.grouping'=gene.grouping, 'sampleset.grouping'=sampleset.grouping,
                'plot.title'=plot.title, 'hide.legend'=hide.legend, 'sort.method' = sort.method,
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

# PlotMaster H.A.K. [Hugs And Kisses] Felix Hu, Justin Chen ---------------------------------------------

# Returns a palette of size n
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

# Returns the first or second token of a split string by ';'
labSplit<-function(str, front = TRUE, delimeter = ";"){
  t<-unlist(strsplit(str, delimeter))
  if(front) {
    return(t[1])
  } else {
    return(t[2])
  }
}

# Return generated plots to act as side bars
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

# Returns indices to sort data by method specified in the argument.
# Sorts dat,  xy = 1 sorts x axis, xy = 2 sorts y axis
indices <- function(dat, method = "input", xy = 1){

  dataCl = dat  # Makes a copy of dat                                   
  
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

# Returns a gtable_frame object used for panel alignment from Baptiste's gridExtra library.
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

# Return a dendrogram ggPlot
# xy = 1 operates on x axis, xy = 2 operates on y axis
create.dendrogram <- function(dat, xy = 1) {
  if(xy == 1) {
   dd <- dendro_data(as.dendrogram(hclust(dist(t(dat)))))
   dendrogram.plot <- ggplot() + geom_segment(data = dd$segments, aes(x=x, y=y, xend=xend, yend=yend), color="gray50") + scale_y_reverse() + theme_dendro() + theme(plot.margin=unit(c(0,2,0,1), "npc"))
  } else {
    dd <- dendro_data(as.dendrogram(hclust(dist(dat))))
    dendrogram.plot <- ggplot(segment(dd)) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend), color="gray50") +  coord_flip() + scale_y_reverse() + theme_dendro() + theme(plot.margin=unit(c(0,0,0,0), "npc"))
  }
  return(dendrogram.plot)
}

met = args$sort.method # Method of sorting data

# Reformat 'data' samples
data <- data[,indices(data, method = met)]

# Reformat 'data' genes
data <- data [indices(data, xy = 2, method = met),]

# Read plot data
plotData<-as.data.frame(melt(as.matrix(data)))
names(plotData) = c('Gene', 'Sample', 'Variant')
#Keep track of original indices of plotData
plotData$ind <- c(1:nrow(plotData))

# Gene grouping
if(!is.na(args$gene.grouping)){
  names(gene_categories) = c("Gene", "Gene_Category")
  
    # Error checking
    if(length(gene_categories$Gene) < length(unique(plotData$Gene))){
      stop("Your plotting data and grouping manual do not match: \"", c(setdiff(as.character(gene_categories$Gene), as.character(unique(plotData$Gene))), setdiff(as.character(unique(plotData$Gene)), as.character(gene_categories$Gene))), "\" fail to overlap")
    }
  
  plotData<-merge(x = plotData, y = gene_categories, all.x = TRUE)
}

# Sampleset grouping
if(!is.na(args$sampleset.grouping)){
  names(sample_categories) = c("Sample", "Sample_Category")
  
    # Error checking
    if(length(sample_categories$Sample) < length(unique(plotData$Sample))){
      stop("Your plotting data and sampleset grouping manual do not match: \"", c(setdiff(as.character(sample_categories$Sample), as.character(unique(plotData$Sample))), setdiff(as.character(unique(plotData$Sample)), as.character(sample_categories$Sample))), "\" fail to overlap")
    }  
  
  plotData<-merge(x = plotData, y = sample_categories, all.x = TRUE)
}

# Sort plotData according to original indices
plotData <- plotData[order(plotData$ind),]

# Prepare y-axis labels if gene grouping is true
if(!is.na(args$gene.grouping)){
  labY = unique(sort(interaction(plotData$Gene_Category, plotData$Gene, lex.order = FALSE, sep = ";")))
  lY2<-unlist(lapply(labY, function(x) labSplit(as.character(x))))
  lY<-unlist(lapply(labY, function(x) labSplit(as.character(x), FALSE)))
}

# Prepare x-axis labels if sampleset grouping is true
if(!is.na(args$sampleset.grouping)){
  labX = unique(sort(interaction(plotData$Sample_Category, plotData$Sample, lex.order = FALSE, sep = ";")))
  lX2<-unlist(lapply(labX, function(x) labSplit(as.character(x))))
  lX<-unlist(lapply(labX, function(x) labSplit(as.character(x),front = FALSE)))
}

# Heatmap Plot
# Tell ggPlot that plotData$Sample is a factor
plotData$Sample <- as.character(plotData$Sample)
plotData$Sample <- factor(plotData$Sample, levels = unique(plotData$Sample))

# Tell ggPlot that plotData$Gene is a factor
plotData$Gene <- as.character(plotData$Gene)
plotData$Gene <- factor(plotData$Gene, levels = unique(plotData$Gene))

plot <- ggplot(plotData, aes(x=Sample, y=Gene))

if(args$gradiant.scale){
  plot = plot + geom_tile(aes(fill = Variant)) + scale_fill_gradient(low = "white", high = "black", name = "Key") + theme_bw()
} else {
  plot = plot + geom_tile(aes(fill = as.character((Variant)))) + scale_fill_manual(values = c("white", getPalette(length(unique(plotData$Variant)))), name = "Key") + theme_bw()
}

# Format x axis labels
if(args$skip.x.label == TRUE){
  textx = element_blank()
} else {
  textx = element_text(angle=-90, hjust=0, vjust=1, color = "black")
}

# Format y axis labels
if(args$skip.y.label == TRUE){
  texty = element_blank()
} else {
  texty = element_text(color = "black")
}

plot = plot + theme(axis.title.x = element_blank(), axis.text.x = textx, axis.title.y = element_blank(), axis.text.y = texty, panel.grid = element_blank(),
        axis.ticks = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        plot.margin = unit(c(0.01, 0.01, 0, 0),"npc"), legend.direction = "horizontal", legend.position = "right")

# Legend
if(args$hide.legend == TRUE){
  plot = plot + guides(fill = FALSE)
}
  
# Create top bar for sampleset grouping
if(!is.na(args$sampleset.grouping)){
  # Get top bar coordinates
  indX<-c(which(diff(as.numeric(as.factor(lX2)))!=0))
  indX = c(indX, length(lX2))
}

# Create side bar for gene grouping
if(!is.na(args$gene.grouping)){
  # Get side bar coordinates
  indY<-c(which(diff(as.numeric(as.factor(lY2)))!=0))
  indY = c(indY, length(lY2))
}

# Build heatmap plot
plot.gt.noguides <- ggplot_gtable(ggplot_build(plot + guides(fill = FALSE)))

#plot.gt <- gtable_add_rows(plot.gt, unit(5,"cm"))
#plot.gt <- gtable_add_grob(plot.gt, ggplotGrob(create.dendrogram(data.original, xy = 1)), t = nrow(plot.gt), l=4, b=nrow(plot.gt), r=4)
#plot.gt <- gtable_add_cols(plot.gt, unit(5, "cm"), pos = 0)

plot.gt.guides <- ggplot_gtable(ggplot_build(plot))
plot.gt.noguides$layout$clip[plot.gt.noguides$layout$name == "panel"] <- "off"
legend.plot <- gtable_filter(plot.gt.guides, "guide-box")

# Build topbar plot
topbar.gp <- bars(dat = plotData, RColBrewPalette = "Set2", isTop = T)
topbar.gp.noguides = topbar.gp + guides(fill = FALSE)
topbar.gt.noguides <- ggplot_gtable(ggplot_build(topbar.gp.noguides))
topbar.gt <- ggplot_gtable(ggplot_build(topbar.gp))
topbar.gt.noguides$layout$clip[topbar.gt.noguides$layout$name == "panel"] <- "off"
legend.topbar <- gtable_filter(topbar.gt, "guide-box")
topbar.gt.noguides$widths = plot.gt.noguides$widths

# Build sidebar plot
sidebar.gp <- bars(dat = plotData, isTop = F) 
sidebar.gp.noguides = sidebar.gp + guides(fill = FALSE)
sidebar.gt <- ggplot_gtable(ggplot_build(sidebar.gp))
sidebar.gt.noguides <- ggplot_gtable(ggplot_build(sidebar.gp.noguides))
sidebar.gt.noguides$layout$clip[sidebar.gt.noguides$layout$name == "panel"] <- "off"
legend.sidebar <- gtable_filter(sidebar.gt, "guide-box")
sidebar.gt.noguides$heights = plot.gt.noguides$heights

# Create Dendrograms
dendro.bottom.gb = ggplot_build(create.dendrogram(dat = data.original, xy = 1))
dendro.bottom.gb$layout$panel_ranges[[1]]$x.range = ggplot_build(plot)$layout$panel_ranges[[1]]$x.range
dendro.bottom.gt = ggplot_gtable(dendro.bottom.gb)
dendro.bottom.gt$widths = plot.gt.noguides$widths

# Create blank.gt ggPlots
blank.gt <- ggplot_gtable(ggplot_build(ggplot() + geom_blank()))

blank2 <- ggplot_gtable(ggplot_build(ggplot() + geom_blank()))
blank2$heights = dendro.bottom.gt$heights
blank2$widths = dendro.bottom.gt$widths

# Combine legends
legendGrob <- gtable:::cbind_gtable(legend.plot, legend.topbar, "last")
legendGrob2 <- gtable:::cbind_gtable(legendGrob, legend.sidebar, "last")

topbar.gf.noguides <- gtable_frame(topbar.gt.noguides, width = unit(31,"null"), height = unit(1,"null"))
sidebar.gf.noguides <- gtable_frame(sidebar.gt.noguides, width = unit(1, "null"), height = unit(31, "null"))
plot.gf.noguides <- gtable_frame(plot.gt.noguides, width = unit(31, "null"), height = unit(31, "null"))

blank.gf <- gtable_frame(blank.gt, width = unit(10, "null"), height = unit(1, "null"))

plot.topbar.gf <- gtable_frame(rbind(topbar.gf.noguides, plot.gf.noguides), width=unit(31,"null"), height=unit(32,"null"))

blank.sidebar.gf <- gtable_frame(rbind(blank.gf, sidebar.gf.noguides), width = unit(1,"null"), height = unit(32,"null"))

graph.bars <- gtable_frame(cbind(plot.topbar.gf, blank.sidebar.gf), width = unit(32, "null"), height = unit(32, "null"))

# Combine graph.bars with bottom dendrogram
#blank2.fg <- gtable_frame(blank.gt, width = unit(1, "null"), height = unit(5, "null"))

#------------------------------------------------------------------------------------------------------------
# Save
cat(sprintf("Saved to %s\n", args$plot.fn))
pdf(args$plot.fn, height=args$height, width=args$width)

#grid::grid.draw(create.dendrogram(data, xy = 1))
grid::grid.draw(graph.bars)
#grid.arrange(graph.bars, bottom.dendrogram.gt, side.dendrogram.gt, widths = c(2, 31, 1), heights = c(1, 31, 2), layout_matrix = rbind(c(NA,1,1),
#                                                                                                                                        c(3,1,1),
#                                                                                                                                        c(NA,2,NA)))

#Closer answers
#main.plot.gb <- arrangeGrob(topbar.gt.noguides, blank.gt, plot.gt, sidebar.gt.noguides, widths = c(45, 1), heights = c(1, 20), ncol = 2, nrow = 2) #layout_matrix = rbind(c(1, NA), c(1, 2)))

#main.plot.gb <- arrangeGrob(topbar.gt.noguides, plot.gt, sidebar.gt.noguides, legendGrob2, widths = c(45, 1), heights = c(1, 20, 1), layout_matrix = rbind(c(1, NA),
#                                                                                                                    c(2, 3),
#                                                                                                                    c(4, 4)))



#Closest answer
#grid.arrange(graph.bars, legendGrob2, heights = c(8, 1), nrow = 2, newpage = FALSE)
  
grid.newpage()
grid.draw(dendro.bottom.gt)

grid.newpage()
grid.draw(rbind(plot.gf.noguides, gtable_frame(dendro.bottom.gt)))

grid.newpage()
grid.arrange(plot.gf.noguides, dendro.bottom.gt, nrow = 2, heights = c(2,1))

dev.off()

