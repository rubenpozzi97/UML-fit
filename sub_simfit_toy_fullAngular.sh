#!/bin/bash

# Create directory for log files
if [ ! -d logs_parSub ]; then mkdir logs_parSub; fi

# Creation of the submit HTCondor file
cat << EOF > temp_sub_simfit_toy_fullAngular_oneBin.sub
Executable  = run_simfit_toy_fullAngular.sh
seed        = ( \$(ProcId) / 6 ) + 1
bin         = \$(ProcId) % 6
Arguments   = \$INT(seed) \$INT(bin)
Log         = logs_parSub/sub_\$(ClusterId).log
Output      = logs_parSub/simfit_toy_fullAngular_randLikelihood_\$INT(seed)_\$INT(bin).out
Error       = logs_parSub/simfit_toy_fullAngular_randLikelihood_\$INT(seed)_\$INT(bin).err
+JobFlavour = "testmatch"
EOF

if [ "${USER}" == "fiorendi" ]; then
    echo '+AccountingGroup = "group_u_CMST3.all"'>>temp_sub_simfit_toy_fullAngular_oneBin.sub
fi
echo 'Queue 600'>>temp_sub_simfit_toy_fullAngular_oneBin.sub

# Compilation, submission and file removal
if make simfit_toy_fullAngular
then condor_submit temp_sub_simfit_toy_fullAngular_oneBin.sub
fi
rm temp_sub_simfit_toy_fullAngular_oneBin.sub
