#!/bin/bash

oldPressure=$(awk -F_ 'NR==2 {print $2}' run_python_tekmc)
echo "Previous run of pressure: ${oldPressure}"

for mixture in CO2
do
    gasA=$(echo ${mixture} | cut -d "_" -f 1)
    #gasB=$(echo ${mixture} | cut -d "_" -f 3)

    gasList="${gasA} ${gasB}"

    for gas in ${gasList}
    do
        for mer in 50
        do
            for pressure in 2e6
            do
                    save_dir=production_5us
                    mkdir -p ${save_dir}/logs

                    spacings=$(awk -F "[][{}]" 'BEGIN {OFS=","} /minimum MSE/ {print $2-0.01, $2, $2+0.01}' optimization/${mer}mer_s_${mixture}_308K_${pressure}Pa_${gas}/normalized_MSE.dat)

                    echo "Production spacings for ${mer}mer with ${gas} at ${pressure} Pa: ${spacings}" | tee -a ${save_dir}/logs/${mer}mer_${mixture}_${pressure}_${gas}.log

                    sed -i "7c topology_file       = 'system/${mer}mer_s_${mixture}_308K_${pressure}Pa_${gas}.gro'" inputs.py
                    sed -i "8c trajectory_file     = 'system/${mer}mer_s_${mixture}_308K_${pressure}Pa_${gas}.dcd'" inputs.py
                    sed -i "9c atom_name           = '${gas}'" inputs.py

                    sed -i "11c spacings_list       = [${spacings}]" inputs.py
                    sed -i "12c n_walks             = 500" inputs.py

                    sed -i "22c save_dir            = '${save_dir}'" inputs.py
                    sed -i "23c n_steps             = 100" inputs.py

                    sed -i "2c #PBS -N ${mer}mer_${gas}_${pressure}" run_python_tekmc
                    sed -i "10c python3 run.py >> ${save_dir}/logs/${mer}mer_${mixture}_${pressure}_${gas}.log 2>&1" run_python_tekmc

                    sleep 1s
                    qsub run_python_tekmc
                    sleep 1s

                    oldPressure="${pressure}"
            done
        done
    done
done
