#!/bin/bash
mesh_name="$1"
set --
source /usr/lib/openfoam/openfoam2412/etc/bashrc
set -euo pipefail
case "$mesh_name" in B00_R07|C01_R07|C02_R07) ;; *) exit 2;; esac
cd "/workspace/scratch/tbm-mesh-repair-20261125/$mesh_name"
date -u +%FT%TZ > repair-check-start-UTC.txt
checkMesh -allTopology > log.checkMesh 2>&1
checkMesh -allTopology -allGeometry > log.checkMesh-additional 2>&1 || true
date -u +%FT%TZ > repair-check-end-UTC.txt
