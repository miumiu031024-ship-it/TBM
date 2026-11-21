#!/bin/bash
ramp_case_name=$1
ramp_start=$2
ramp_end=$3
set --
source /usr/lib/openfoam/openfoam2412/etc/bashrc
set -euo pipefail
cd /workspace/scratch/tbm-steady-ramp-20261121/"$ramp_case_name"
test ! -e "log.simpleFoam-$ramp_start-$ramp_end"
if test ! -d processor0; then
 if test "$ramp_case_name" = F02_second_order_restore; then
  # Preserve the exact four-way mesh partition and local restart values.
  for ramp_rank in 0 1 2 3; do
   mkdir "processor$ramp_rank"
   cp -a "../F01_upwind_long/processor$ramp_rank/constant" "processor$ramp_rank/constant"
   cp -a "../F01_upwind_long/processor$ramp_rank/$ramp_start" "processor$ramp_rank/$ramp_start"
  done
  echo 'Preserved F01 four-way decomposition and processor restart fields' > log.decomposePar
 else
  decomposePar -time "$ramp_start" > log.decomposePar 2>&1
 fi
fi
foamDictionary system/controlDict -entry endTime -set "$ramp_end" > "log.control-$ramp_end" 2>&1
mpirun --mca btl_vader_single_copy_mechanism none -np 4 simpleFoam -parallel > "log.simpleFoam-$ramp_start-$ramp_end" 2>&1
reconstructPar -time "$((ramp_end-100)),$ramp_end" -fields '(U p phi k omega nut)' > "log.reconstruct-$ramp_end" 2>&1
