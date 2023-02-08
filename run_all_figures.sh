
FIGURE_FOLDERS=$(ls figures)
CURRENT_DIR=${PWD}

for fig in ${FIGURE_FOLDERS}
do
	cd figures/${fig}
	echo ${PWD}

	run_script="run_${fig}.sh"
	./${run_script}

	cd ${CURRENT_DIR}
done