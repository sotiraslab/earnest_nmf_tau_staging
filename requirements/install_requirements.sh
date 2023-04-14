#!/bin/bash

echo ""
echo "install R requirements"
echo '----------------------'
SCRIPT="install_r_requirements.R"
Rscript -e "source('${SCRIPT}', echo=T)"

echo ""
echo "install python requirements"
echo '---------------------------'
pip install -r requirements_python.txt