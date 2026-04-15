#!/bin/bash -x

echo "Loading cluster environment." 
source ~/MinimalFunctionCache/src/minimalFunctionCache/parametericCutGen/cluster_enviroment.sh
echo "Setting up file system."

if [ ! -d $EXP_TEMP ]; then 
mkdir $EXP_TEMP
fi

if [ ! -d $MODEL_FILES ]; then 
   mkdir $MODEL_FILES
fi

cleanup_file_sys() {
  rm -rf $EXP_TEMP
  rm -rf $MODEL_FILES
  rm -rf $EXP_FILES_TO_RUN
}

if [ ! -d $DATA_OUT_FILE_PATH ]; then 
   mkdir $DATA_OUT_FILE_PATH
else
    if [ (ls -A $DATA_OUT_FILE_PATH) ]; then 
         echo "Data out path contains previously written data. This data will be over written."
         echo "Are you sure you wish to continue [y/n]?"
         read cont
         valid_input=0
         while [ $valid_input -eq 0 ]
         do
             if [ "$cont" ='y' ]; then
             	echo "Removing old data."
             	rm -rf $DATA_OUT_FILE_PATH
             	valid_input=1
             	break
             elif [ $cont=='n' ]; then 
             	echo "Aborting. Cleaning up."
             	cleanup_file_sys()
             	exit 0
             else
                echo "Please input a valid input [y/n]"
                read cont
             fi
         done
fi

if [ ! -d $TRIAL_FILES_TO_RUN ]; then
  mkdir $TRIAL_FILES_TO_RUN
else
   echo "Warning! Previous files to run folder found. Removing."
   rm -rf $TRIAL_FILES_TO_RUN
fi

echo "Checking for apptainer."

module load apptainer
if [ ! -f optimal_cut.sif ]; then
  echo ""
  apptainer build optimal_cut.sif $APPTAINER_DEF_PAT
fi

echo "Configuring trials."
# exp_setup.sh calls # write_trials() with path information following the apptainer.
# write_trials() generates files; trials_env.sh, run_trial.sh, trial_$NUMBER.sh
sbatch --partition=$PARTITION --account=$CLUSTER_ACCOUNT --ntasks=1 --cpus-per-task=1 --time=$SETUP_TIME:00 exp_setup.sh

source "$EXP_TEMP/trials_env.sh"

echo "Running trials."
sbatch --array=0-$NUM_TRIALS--partition=$PARTITION --account=$CLUSTER_ACCOUNT --ntasks=1 --cpus-per-task=1 --time=$TRIAL_RUN_TIME:00 run_trial.sh
