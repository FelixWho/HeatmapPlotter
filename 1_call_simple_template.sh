DAT="testchrom.short1000.txt"
OUT=“output.pdf"

Rscript HeatmapPlotter.R -f -G -P -x -y -m clust -s -v -w 15 -h 15 -t TestPlot $DAT $OUT


