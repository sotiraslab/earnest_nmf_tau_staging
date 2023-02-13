#!/bin/bash

echo ""
echo "----------"
echo "DESCRIPTIVE STATISTOCS"
echo "----------"

echo ""
echo "Running descriptive stats..."
Rscript -e "source('descriptive_stats.R', echo=T)"
echo "Done."

