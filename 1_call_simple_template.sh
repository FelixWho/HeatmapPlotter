DAT="testchrom.short1000.txt"
OUT="cluster.pdf"

Rscript HeatmapPlotter.R -x -y -f -m clust -s  -v -w 15 -h 15  -t TestPlot.Dendro $DAT $OUT

