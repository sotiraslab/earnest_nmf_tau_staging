#!/bin/bash

echo ""
echo "----------"
echo "FIGURE 6"
echo "----------"

echo ""
echo "WARNING: This must be run after table2 has been run!"
echo ""

echo ""
echo "Running plots of gene expression..."
SCRIPT="plot_expression_wightman.R"
Rscript -e "source('${SCRIPT}', echo=T)"
echo "Done."