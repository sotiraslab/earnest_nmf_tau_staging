
FIGURE_FOLDERS=$(ls supplement)
CURRENT_DIR=${PWD}

echo ""
echo "NOTE:"
echo "  Some supplemental materials are produced and saved"
echo "  in the nmf_tau/figures folder!"
echo ""

for fig in ${FIGURE_FOLDERS}
do
	cd supplement/${fig}

	run_script="run.sh"
	if [[ -e ${run_script} ]]
	then
		./${run_script}
	fi
	
	cd ${CURRENT_DIR}
done