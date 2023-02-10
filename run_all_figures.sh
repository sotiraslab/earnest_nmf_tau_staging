
FIGURE_FOLDERS=$(ls figures)
CURRENT_DIR=${PWD}

for fig in ${FIGURE_FOLDERS}
do
	cd figures/${fig}

	run_script="run_${fig}.sh"
	if [[ -e ${run_script} ]]
	then
		./${run_script}
	fi
	
	cd ${CURRENT_DIR}
done