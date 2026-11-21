#!/bin/bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
set -euo pipefail
for name in S00_reference S10_pressure1 S20_velocity015 S30_upwind; do
 cd /workspace/scratch/tbm-cause-20261121/"$name"
 test ! -e log.simpleFoam
 decomposePar -force -time 500 > log.decomposePar 2>&1
 mpirun --mca btl_vader_single_copy_mechanism none -np 4 simpleFoam -parallel > log.simpleFoam 2>&1
 reconstructPar -time '550,600' -fields '(U p k omega nut phi)' > log.reconstructPar 2>&1
 foamToVTK -time 600 -fields '(U p)' > log.foamToVTK 2>&1
 echo "$name complete"
done
