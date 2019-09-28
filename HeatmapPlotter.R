library(ggdendro)
library(plyr)
library(tidyr)
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
    side.dendrogram = get_bool_arg(args, "-G") 
    bottom.dendrogram = get_bool_arg(args, "-P") 
    plot.title = get_val_arg(args, "-t", "Plot")
    hide.legend = get_bool_arg(args, "-l") #TODO: FIX
    gradiant.scale = get_bool_arg(args, "-s")
    sort.method = as.character(get_val_arg(args, "-m", "alphabetgrouping"))
    data.input.longFormat = get_bool_arg(args, "-f")
    
    # mandatory positional arguments.  These are popped off the back of the array, last one listed first.
    plot.fn = args[length(args)];               args = args[-length(args)]
    data.fn = args[length(args)]; args = args[-length(args)]
    
    
    val = list( 'plot.fn'=plot.fn, 'verbose'=verbose, 'data.fn'=data.fn,
                'height'=height, 'width'=width, 'skip.x.label'=skip.x.label, 'skip.y.label'=skip.y.label,
                'gene.grouping'=gene.grouping, 'sampleset.grouping'=sampleset.grouping, 'side.dendrogram' = side.dendrogram, 'bottom.dendrogram' = bottom.dendrogram,
                'plot.title'=plot.title, 'hide.legend'=hide.legend, 'sort.method' = sort.method,
                'gradiant.scale'=gradiant.scale, 'sort.method'=sort.method,'data.input.longFormat'=data.input.longFormat)
    
    if (val$verbose) { print(val) }
    return (val)
}

options("width"=180) # useful for debugging
#options("width"=270) 
args = parse_args()

# Read input file(s)
if(args$data.input.longFormat){
  data.long<-read.table(args$data.fn, header = TRUE, sep = "\t", check.names = FALSE, comment.char = "#")
  data <- dcast(data.long, formula = data.long[[1]] ~ data.long[[2]], value.var = colnames(data.long)[3])
  data[[1]] = NULL
  Y.Column = as.matrix(unique(data.long[[2]]))
  X.Column = as.matrix(unique(data.long[[1]]))
  colnames(data) = X.Column
  rownames(data) = Y.Column
} else {
data<-read.table(args$data.fn, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, comment.char = "#") #Reads TSV
}
data.original = data

if(!is.na(args$gene.grouping)){
  gene_categories <- read.table(as.character(args$gene.grouping), header = FALSE, sep = "\t", comment.char = "#")
}
if(!is.na(args$sampleset.grouping)){
  sample_categories <- read.table(as.character(args$sampleset.grouping), header=FALSE, sep = "\t", comment.char = "#")
}

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
    f <- ggplot(val, aes(val, x = Sample, y = Gene)) + geom_tile(aes(fill = as.character(Sample_Category))) + scale_fill_manual(values = c(getPalette(length(unique(val$Sample_Category)))), name = "Gene Key")
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
   dendrogram.plot <- ggplot() + geom_segment(data = dd$segments, aes(x=x, y=y, xend=xend, yend=yend), color="gray50") + scale_y_reverse() + theme_dendro() + theme(plot.margin=unit(c(0,0,0,0), "npc"))
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
data <- data[indices(data, xy = 2, method = met),]

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
        axis.ticks = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.5,),
        plot.margin = unit(c(0, 0, 0, 0),"npc"), legend.direction = "horizontal", legend.position = "right")

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

# Build heatmap
plot.gt.noguides <- ggplot_gtable(ggplot_build(plot + guides(fill = FALSE))) # plot without legend
plot.gt.noguides$layout$clip[plot.gt.noguides$layout$name == "panel"] <- "off"
plot.gf.final <- gtable_frame(plot.gt.noguides, width = unit(31, "null"), height = unit(31, "null"))

if(!args$hide.legend){
  plot.gt.guides <- ggplot_gtable(ggplot_build(plot)) # plot with legend
  legend.grob <- gtable_filter(plot.gt.guides, "guide-box") # rip the legend off plot.gt.guides
}

if(!is.na(args$sampleset.grouping)) {
  # Build topbar plot
  topbar.gp <- bars(dat = plotData, RColBrewPalette = "Set2", isTop = T)
  topbar.gp.noguides = topbar.gp + guides(fill = FALSE)
  topbar.gt.noguides <- ggplot_gtable(ggplot_build(topbar.gp.noguides))
  topbar.gt.noguides$widths = plot.gt.noguides$widths
  topbar.gt.noguides$layout$clip[topbar.gt.noguides$layout$name == "panel"] <- "off"
  topbar.gf.noguides <- gtable_frame(topbar.gt.noguides, width = unit(31,"null"), height = unit(1,"null"))
  
  plot.gf.final = gtable_frame(rbind(topbar.gf.noguides, plot.gf.final), width=unit(31,"null"), height=unit(32,"null"))
  
  if(!args$hide.legend){
    topbar.gt <- ggplot_gtable(ggplot_build(topbar.gp))
    legend.topbar <- gtable_filter(topbar.gt, "guide-box") # rip the legend off of topbar.gt
    legend.grob = gtable:::cbind_gtable(legend.grob, legend.topbar, "last")
  }
}

if(!is.na(args$gene.grouping)) {
  # Build sidebar plot
  sidebar.gp <- bars(dat = plotData, isTop = F) 
  sidebar.gp.noguides = sidebar.gp + guides(fill = FALSE)
  sidebar.gt.noguides <- ggplot_gtable(ggplot_build(sidebar.gp.noguides))
  sidebar.gt.noguides$heights = plot.gt.noguides$heights
  sidebar.gt.noguides$layout$clip[sidebar.gt.noguides$layout$name == "panel"] <- "off"
  sidebar.gf.noguides <- gtable_frame(sidebar.gt.noguides, width = unit(1, "null"), height = unit(31, "null"))
  
  if(!is.na(args$sampleset.grouping)){
    blank.gt <- ggplot_gtable(ggplot_build(ggplot() + geom_blank()))
    blank.gt$heights = topbar.gt.noguides$heights
    blank.gt$widths = topbar.gt.noguides$widths
    blank.gf = gtable_frame(blank.gt, width = unit(1, "null"), height = unit(1, "null"))
    
    sidebar.gf <- gtable_frame(rbind(blank.gf, sidebar.gf.noguides), width = unit(1,"null"), height = unit(32,"null"))
    plot.gf.final = gtable_frame(cbind(plot.gf.final, sidebar.gf), width = unit(32, "null"), height = unit(32, "null"))
  
  } else {
    sidebar.gf <- sidebar.gf.noguides
    plot.gf.final = gtable_frame(cbind(plot.gf.final, sidebar.gf), width = unit(32, "null"), height = unit(31, "null"))
  }

  if(!args$hide.legend){
    sidebar.gt <- ggplot_gtable(ggplot_build(sidebar.gp))
    legend.sidebar <- gtable_filter(sidebar.gt, "guide-box") # rip the legend off of sidebar.gt
    legend.grob = gtable:::cbind_gtable(legend.grob, legend.sidebar, "last")
  }
}

if(args$side.dendrogram){
  # Build side dendrogram
  dendro.side.gb = ggplot_build(create.dendrogram(dat = data.original, xy = 2)+theme(panel.background = element_rect(fill = 'red', colour = 'red')))
  dendro.side.gb$layout$panel_ranges[[1]]$y.range = ggplot_build(plot)$layout$panel_ranges[[1]]$y.range
  dendro.side.gt = ggplot_gtable(dendro.side.gb)
  dendro.side.gt$heights = plot.gt.noguides$heights
  dendro.side.gf = gtable_frame(dendro.side.gt, width = unit(1, "null"), height = unit(31, "null"))
  
  if(!is.na(args$sampleset.grouping)){
  blank.gt <- ggplot_gtable(ggplot_build(ggplot() + geom_blank() +theme(plot.margin=unit(c(0,0,0,0), "npc"),panel.background = element_rect(fill = 'green', colour = 'red'))))
  blank.gt$heights = dendro.side.gt$heights
  blank.gt$widths = dendro.side.gt$widths
  blank.gf = gtable_frame(blank.gt, width = unit(1, "null"), height = unit(1, "null"))
  
  dendro.side.final <- gtable_frame(rbind(blank.gf, dendro.side.gf), width = unit(1,"null"), height = unit(32,"null"))
    if(!is.na(args$gene.grouping)){
    plot.gf.final = gtable_frame(cbind(dendro.side.final, plot.gf.final), width = unit(33,"null"), height = unit(32,"null"))
    } else {
      plot.gf.final = gtable_frame(cbind(dendro.side.final, plot.gf.final), width = unit(32,"null"), height = unit(32,"null"))
    }
  } else if(!is.na(args$gene.grouping)){
      plot.gf.final = gtable_frame(cbind(dendro.side.gf, plot.gf.final), width = unit(33,"null"), height = unit(31,"null"))
    } else {
      plot.gf.final = gtable_frame(cbind(dendro.side.gf, plot.gf.final), width = unit(32,"null"), height = unit(31,"null"))
    }
  
}

if(args$bottom.dendrogram){
  # Build bottom dendrogram
  dendro.bottom.gb = ggplot_build(create.dendrogram(dat = data.original, xy = 1) +theme(panel.background = element_rect(fill = 'red', colour = 'red')))
  dendro.bottom.gb$layout$panel_ranges[[1]]$x.range = ggplot_build(plot + guides(fill = FALSE))$layout$panel_ranges[[1]]$x.range
  dendro.bottom.gt = ggplot_gtable(dendro.bottom.gb)
  dendro.bottom.gt$widths = plot.gt.noguides$widths
  dendro.bottom.final = gtable_frame(dendro.bottom.gt, width = unit(31, "null"), height = unit(1, "null"))
  
  if(!is.na(args$gene.grouping)){
    blank.gt <- ggplot_gtable(ggplot_build(ggplot() + geom_blank()+ theme(plot.margin=unit(c(0,0,0,0), "npc"),panel.background = element_rect(fill = 'green', colour = 'red'))))
    blank.gt$heights = dendro.bottom.gt$heights
    blank.gt$widths = dendro.bottom.gt$widths
    blank.gf = gtable_frame(blank.gt, width = unit(1, "null"), height = unit(1, "null"))
    
    dendro.bottom.final = gtable_frame(cbind(dendro.bottom.final, blank.gf), width = unit(32,"null"), height = unit(1,"null"))
    if(args$side.dendrogram){
      dendro.bottom.final = gtable_frame(cbind(blank.gf, dendro.bottom.final), width = unit(33,"null"), height = unit(1,"null"))
    }
  } else if(args$side.dendrogram){
    blank.gt <- ggplot_gtable(ggplot_build(ggplot() + geom_blank()+ theme(plot.margin=unit(c(0,0,0,0), "npc"),panel.background = element_rect(fill = 'green', colour = 'red'))))
    blank.gt$heights = dendro.bottom.gt$heights
    blank.gt$widths = dendro.bottom.gt$widths
    blank.gf = gtable_frame(blank.gt, width = unit(1, "null"), height = unit(1, "null"))
    
    dendro.bottom.final = gtable_frame(cbind(blank.gf, dendro.bottom.final), width = unit(32,"null"), height = unit(1,"null"))
  }
  
  if(!is.na(args$sampleset.grouping)){
    if(args$side.dendrogram & !is.na(args$gene.grouping)) {
      plot.gf.final = gtable_frame(rbind(plot.gf.final, dendro.bottom.final), width = unit(33,"null"), height = unit(33,"null"))
    } else if(args$side.dendrogram | !is.na(args$gene.grouping)){
      plot.gf.final = gtable_frame(rbind(plot.gf.final, dendro.bottom.final), width = unit(32,"null"), height = unit(33,"null"))
    }
  } else {
    if(args$side.dendrogram & !is.na(args$gene.grouping)) {
      plot.gf.final = gtable_frame(rbind(plot.gf.final, dendro.bottom.final), width = unit(33,"null"), height = unit(32,"null"))
    } else if(args$side.dendrogram | !is.na(args$gene.grouping)){
      plot.gf.final = gtable_frame(rbind(plot.gf.final, dendro.bottom.final), width = unit(32,"null"), height = unit(32,"null"))
    }
  }
}

# Title
if (!is.null(args$plot.title)) {
  title.ggp = textGrob(args$plot.title)
  plot.gf.final = arrangeGrob(title.ggp, plot.gf.final, heights = c(0.025, 0.975), ncol=1, nrow=2)
}

#------------------------------------------------------------------------------------------------------------
# Save
cat(sprintf("Saved to %s\n", args$plot.fn))
pdf(args$plot.fn, height=args$height, width=args$width)
                                                                                                                              
grid.draw(plot.gf.final)


dev.off()

