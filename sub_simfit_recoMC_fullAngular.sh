#!/bin/bash

# Create directory for log files
if [ ! -d logs_parSub ]; then mkdir logs_parSub; fi

make simfit_recoMC_fullAngular_parSub

# Creation of the submit HTCondor file
cat << EOF > temp_sub_simfit_recoMC_fullAngular_oneBin.sub
Executable  = run_simfit_recoMC_fullAngular.sh
bin         = 0
f1          = 100
f4          = 10
f5          = 10
b1          = ( \$(ProcId) / 7 ) * 4 + 36
b4          = ( \$(ProcId) % 7 ) * 4 + 36
b5          = 1
m1          = 0
m4          = 0
m5          = 0
Arguments   = \$INT(f1) \$INT(f4) \$INT(f5) \$INT(b1) \$INT(b4) \$INT(b5) \$INT(m1) \$INT(m4) \$INT(m5) \$INT(bin)
Log         = logs_parSub/sub_\$(ClusterId).log
Output      = logs_parSub/simfit_recoMC_fullAngular_parSub11_\$INT(f1)_\$INT(f4)_\$INT(f5)_\$INT(b1)_\$INT(b4)_\$INT(b5)_\$INT(m1)_\$INT(m4)_\$INT(m5)_\$INT(bin).out
Error       = logs_parSub/simfit_recoMC_fullAngular_parSub11_\$INT(f1)_\$INT(f4)_\$INT(f5)_\$INT(b1)_\$INT(b4)_\$INT(b5)_\$INT(m1)_\$INT(m4)_\$INT(m5)_\$INT(bin).err
+JobFlavour = "testmatch"
EOF

if [ "${USER}" == "fiorendi" ]; then
    echo '+AccountingGroup = "group_u_CMST3.all"'>>temp_sub_simfit_recoMC_fullAngular_oneBin.sub
fi
echo 'Queue 49'>>temp_sub_simfit_recoMC_fullAngular_oneBin.sub

# Submission and file removal
condor_submit temp_sub_simfit_recoMC_fullAngular_oneBin.sub
rm temp_sub_simfit_recoMC_fullAngular_oneBin.sub
