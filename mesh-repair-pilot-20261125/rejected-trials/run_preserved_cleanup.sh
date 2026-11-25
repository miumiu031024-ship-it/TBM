#!/bin/bash
task_mesh_name="$1"
set --
source /usr/lib/openfoam/openfoam2412/etc/bashrc
set -euo pipefail
root=/workspace/scratch/tbm-mesh-repair-20261125
case "$task_mesh_name" in B00_R03|C01_R03|C02_R03) ;; *) exit 2;; esac
cd "$root/$task_mesh_name"
test ! -e log.removeFaces
date -u +%FT%TZ > repair-start-UTC.txt
if test -f selected-repairFaces.input; then
 mv 0 initial-conditions-preserved
 mkdir 0
 set +e
 /workspace/scratch/tbm-mesh-repair-20261125/bin/internalFaceAgglomerate -overwrite repairFaces > log.removeFaces 2>&1
 task_cleanup_status=$?
 set -e
 mv 0 mesh-utility-start-output
 mv initial-conditions-preserved 0
 test "$task_cleanup_status" -eq 0
else
 echo 'No cells below the determinant criterion; identical quality-driven rule selects no repair.' > log.removeFaces
fi
checkMesh -allTopology > log.checkMesh 2>&1
checkMesh -allTopology -allGeometry > log.checkMesh-additional 2>&1 || true
date -u +%FT%TZ > repair-end-UTC.txt
echo "$task_mesh_name local cleanup and inspections finished"
