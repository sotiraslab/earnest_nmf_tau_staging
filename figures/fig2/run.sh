#!/bin/bash

echo ""
echo "----------"
echo "W-Scoring"
echo "----------"

echo ""
echo "Running W-scoring & heatmaps..."
SCRIPT="wscore_heatmaps.R"
Rscript "${SCRIPT}"
echo "Done."

echo ""
echo "Running winner take all image..."
SCRIPT="winner_take_all.R"
Rscript "${SCRIPT}"
echo "Done."

echo ""
echo "Running W-score brain maps..."
python wscore_brainmaps.py
echo "Done."