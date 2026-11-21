#!/bin/bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
set -euo pipefail
for name in T30_dt2e5 T40_dt1e5; do
 cd /workspace/scratch/tbm-cause-20261121/"$name"
 test ! -e log.pimpleFoam
 decomposePar > log.decomposePar 2>&1
 mpirun --mca btl_vader_single_copy_mechanism none -np 4 pimpleFoam -parallel > log.pimpleFoam 2>&1
 reconstructPar -latestTime -fields '(U p phi)' > log.reconstructPar 2>&1
 foamToVTK -latestTime -fields '(U p)' > log.foamToVTK 2>&1
 echo "$name complete"
done
