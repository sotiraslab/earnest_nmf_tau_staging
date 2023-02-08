#!/bin.bash

echo ""
echo "----------"
echo "FIGURE 2"
echo "----------"

echo ""
echo "Running creation of component images..."
Rscript -e "source('create_8ptc_ggseg.R', echo=T)"
echo "Done."