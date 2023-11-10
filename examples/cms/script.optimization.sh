#!/bin/bash

oldPressure=$(awk -F_ 'NR==2 {print $2}' run_python_tekmc)
echo "Previous run of pressure: ${oldPressure}"

mixture="CO2"
for gas in CO2
do
    for mer in 50
    do
        for pressure in 2e6
        do
                save_dir=optimization
                mkdir -p ${save_dir}/logs

                spacings=$(awk -F "[][{}]" 'BEGIN {OFS=","} /minimum MSE/ {print $2-0.03, $2-0.01, $2, $2+0.01, $2+0.03}' initiation/${mer}mer_s_${mixture}_308K_${pressure}Pa_${gas}/normalized_MSE.dat)

                echo "Production spacings for ${mer}mer with ${gas} at ${pressure} Pa: ${spacings}" | tee -a ${save_dir}/logs/${mer}mer_${mixture}_${pressure}_${gas}.log

                sed -i "7c topology_file       = 'system/${mer}mer_s_${mixture}_308K_${pressure}Pa_${gas}.gro'" inputs.py
                sed -i "8c trajectory_file     = 'system/${mer}mer_s_${mixture}_308K_${pressure}Pa_${gas}.dcd'" inputs.py
                sed -i "9c atom_name           = '${gas}'" inputs.py

                sed -i "11c spacings_list       = [${spacings}]" inputs.py
                sed -i "12c n_walks             = 1000" inputs.py

                sed -i "22c save_dir            = '${save_dir}'" inputs.py
                sed -i "23c n_steps             = 10000" inputs.py

                sed -i "2c #PBS -N ${mer}mer_${gas}_${pressure}" run_python_tekmc
                sed -i "10c python3 run.py >> ${save_dir}/logs/${mer}mer_${mixture}_${pressure}_${gas}.log 2>&1" run_python_tekmc

                sleep 1s

                #/home/jeet/installations/anaconda3/bin/python3 run.py >> initiation/logs/${mer}mer_${mixture}_${pressure}_${gas}.log 2>&1
                qsub run_python_tekmc

                sleep 1s

                oldPressure="${pressure}"
        done
    done
done
