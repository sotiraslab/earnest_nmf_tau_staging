#!/bin/bash

echo ""
echo "----------"
echo "CN Demographic table"
echo "----------"

SCRIPT="supplemental_descriptives.R"
Rscript -e "source('${SCRIPT}', echo=T)"