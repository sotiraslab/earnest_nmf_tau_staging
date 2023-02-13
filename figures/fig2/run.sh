#!/bin/bash

echo ""
echo "----------"
echo "W-Scoring"
echo "----------"

echo ""
echo "Running W-scoring & heatmaps..."
SCRIPT="wscore_heatmaps.R"
Rscript -e "source('${SCRIPT}', echo=T)"
echo "Done."

echo ""
echo "Running winner take all image..."
SCRIPT="winner_take_all.R"
Rscript -e "source('${SCRIPT}', echo=T)"
echo "Done."

echo ""
echo "Running W-score brain maps..."
SCRIPT="wscore_brainmaps.R"
Rscript -e "source('${SCRIPT}', echo=T)"
echo "Done."