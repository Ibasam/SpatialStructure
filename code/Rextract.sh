#!/bin/bash

#DIR="results/" # directory where to save results
#cd $DIR


# SIMULATION
echo "BEGINNING OF SIMULATIONS"
echo "PID du processus courant : $$"
STARTTIME=$(date +%s);


scn=(10311111 10341111 10341211 10341222) #example

for s in "${scn[@]}"
do

#STARTTIME=$(date +%s);
Rscript --vanilla DEMOGRAPHY.R $s &
Rscript --vanilla PHENOGENOTYPE.R $s &
#Rscript --vanilla FITNESS.R $s &

echo "Extraction of scenario $s started!" 
#sleep 2 # wait x seconds before starting the next extraction
done # end loop 


#ENDTIME=$(date +%s);
#MINUTES=$(( ($ENDTIME - $STARTTIME) / 60 ));
#echo "Simulations successful! Duration: $MINUTES minutes" 
echo "END OF SIMULATIONS"
