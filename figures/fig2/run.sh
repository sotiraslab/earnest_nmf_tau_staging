#!/bin/bash

echo ""
echo "----------"
echo "PTC VISUALIZATION"
echo "----------"

echo ""
echo "Running creation of component images..."
Rscript -e "source('create_8ptc_ggseg.R', echo=T)"
echo "Done."