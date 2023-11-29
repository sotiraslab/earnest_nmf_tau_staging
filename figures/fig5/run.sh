#!/bin/bash

echo ""
echo "----------"
echo "LOADINGS REGRRESSIONS"
echo "----------"

echo ""
echo "Running Centiloid/loading for ADNI..."
SCRIPT="adni_centiloid_regression.R"
Rscript -e "source('${SCRIPT}', echo=T)"
echo "Done."

echo ""
echo "Running Centiloid/loading for OASIS..."
SCRIPT="oasis_centiloid_regression.R"
Rscript -e "source('${SCRIPT}', echo=T)"
echo "Done."

echo ""
echo "Running volume/loading for ADNI..."
SCRIPT="adni_gm_regression.R"
Rscript -e "source('${SCRIPT}', echo=T)"
echo "Done."

echo ""
echo "Running volume/loading for OASIS..."
SCRIPT="oasis_gm_regression.R"
Rscript -e "source('${SCRIPT}', echo=T)"
echo "Done."