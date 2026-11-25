#!/bin/bash
task_mesh_name="$1"
set --
source /usr/lib/openfoam/openfoam2412/etc/bashrc
set -euo pipefail
root=/workspace/scratch/tbm-mesh-repair-20261125
name="$task_mesh_name"
case "$name" in B00_R01|C01_R01|C02_R01) ;; *) exit 2;; esac
cd "$root/$name"
test ! -e log.snappyHexMesh
date -u +%FT%TZ > mesh-start-UTC.txt
surfaceCheck constant/triSurface/volume.stl > log.surfaceCheck 2>&1
grep -q 'Surface is closed' log.surfaceCheck
blockMesh > log.blockMesh 2>&1
decomposePar -no-fields > log.decomposePar-mesh 2>&1
{ time mpirun --mca btl_vader_single_copy_mechanism none -np 4 snappyHexMesh -parallel -overwrite; } > log.snappyHexMesh 2> log.snappy-time
reconstructParMesh -constant -noZero > log.reconstructParMesh 2>&1
checkMesh -allTopology > log.checkMesh 2>&1
checkMesh -allTopology -allGeometry > log.checkMesh-additional 2>&1 || true
date -u +%FT%TZ > mesh-end-UTC.txt
echo "$name revision mesh and inspections finished"
