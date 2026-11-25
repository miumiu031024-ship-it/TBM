#!/bin/bash
mesh_name="$1"
set --
source /usr/lib/openfoam/openfoam2412/etc/bashrc
set -euo pipefail
case "$mesh_name" in B00_R07|C01_R07|C02_R07) ;; *) exit 2;; esac
cd "/workspace/scratch/tbm-mesh-repair-20261125/$mesh_name"
test ! -e log.simpleFoam20
grep -q 'Mesh OK.' log.checkMesh
if grep -q 'Cells with small determinant' log.checkMesh-additional; then echo 'determinant gate failed'; exit 4; fi
printf '%s\n' "$mesh_name" > ../active-pilot.txt
date -u +%FT%TZ > pilot-start-UTC.txt
cat /sys/fs/cgroup/memory.current > pilot-start-container-memory-bytes.txt
decomposePar -force > log.decomposePar-pilot20 2>&1
date -u +%FT%TZ > solver-start-UTC.txt
TIMEFORMAT='real_seconds %3R user_seconds %3U sys_seconds %3S'
{ time mpirun --mca btl_vader_single_copy_mechanism none -np 4 simpleFoam -parallel > log.simpleFoam20 2>&1; } 2> solver-driver-time.txt
date -u +%FT%TZ > solver-end-UTC.txt
grep -q '^Time = 20' log.simpleFoam20
grep -q '^End' log.simpleFoam20
reconstructPar -time 20 > log.reconstructPar20 2>&1
test -f 20/U || test -f 20/U.gz
test -f 20/p || test -f 20/p.gz
date -u +%FT%TZ > pilot-end-UTC.txt
printf '%s\n' none > ../active-pilot.txt
