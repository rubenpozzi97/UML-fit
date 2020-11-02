#!/bin/bash

nGenMIN=1000000
widthSc=0.05
conf=15

echo ${conf} ${widthSc} ${nGenMIN} >> compConf.list

while read -a line; do

    bin=${line[0]}
    gen1=${line[1]}
    gen2=${line[2]}
    gen3=${line[3]}
    gen4=${line[4]}
    gen5=${line[5]}
    gen6=${line[6]}
    gen7=${line[7]}
    gen8=${line[8]}

    # Create directory for log files
    if [ ! -d logs_parSub/toy_b${bin}_v${conf} ]; then mkdir -p logs_parSub/toy_b${bin}_v${conf}; fi

    # Creation of the submit HTCondor file
    cat << EOF > temp_sub_simfit_toy_fullAngular_oneBin.sub
Executable  = run_simfit_toy_fullAngular.sh
seed        = \$(ProcId) + 1
bin         = ${bin}
Arguments   = \$INT(seed) \$INT(bin) ${gen1} ${gen2} ${gen3} ${gen4} ${gen5} ${gen6} ${gen7} ${gen8} ${nGenMIN} ${widthSc} ${conf}
Log         = logs_parSub/sub_\$(ClusterId).log
Output      = logs_parSub/toy_b${bin}_v${conf}/simfit_${gen1}_${gen2}_${gen3}_${gen4}_${gen5}_${gen6}_${gen7}_${gen8}_${nGenMIN}_${widthSc}_\$INT(seed).out
Error       = logs_parSub/toy_b${bin}_v${conf}/simfit_${gen1}_${gen2}_${gen3}_${gen4}_${gen5}_${gen6}_${gen7}_${gen8}_${nGenMIN}_${widthSc}_\$INT(seed).err
+JobFlavour = "nextweek"
EOF

    if [ "${USER}" == "fiorendi" ]; then
	echo '+AccountingGroup = "group_u_CMST3.all"'>>temp_sub_simfit_toy_fullAngular_oneBin.sub
    fi
    echo 'Queue 50'>>temp_sub_simfit_toy_fullAngular_oneBin.sub
    
    # Compilation, submission and file removal
    if make simfit_toy_fullAngular
    then
	echo "submitting" ${bin} ${gen1} ${gen2} ${gen3} ${gen4} ${gen5} ${gen6} ${gen7} ${gen8}
	condor_submit temp_sub_simfit_toy_fullAngular_oneBin.sub
    fi
    rm temp_sub_simfit_toy_fullAngular_oneBin.sub

done < ../confSF/toy_coord.list
