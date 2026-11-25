#!/bin/bash
set --
source /usr/lib/openfoam/openfoam2412/etc/bashrc
set -euo pipefail
for name in B00 C01 C02; do
  while ! test -s /work/"$name"/case-metadata.json; do sleep 2; done
  cd /work/"$name"
  test ! -e log.snappyHexMesh
  date -u +%FT%TZ > mesh-start-UTC.txt
  surfaceCheck constant/triSurface/volume.stl > log.surfaceCheck 2>&1
  grep -q 'Surface is closed' log.surfaceCheck
  surfaceFeatureExtract > log.surfaceFeatureExtract 2>&1
  blockMesh > log.blockMesh 2>&1
  decomposePar -force -no-fields > log.decomposePar-mesh 2>&1
  { time mpirun --mca btl_vader_single_copy_mechanism none -np 4 snappyHexMesh -parallel -overwrite; } > log.snappyHexMesh 2> log.snappy-resource
  reconstructParMesh -constant -noZero > log.reconstructParMesh 2>&1
  checkMesh -allTopology > log.checkMesh 2>&1
  checkMesh -allTopology -allGeometry > log.checkMesh-additional 2>&1 || true
  date -u +%FT%TZ > mesh-end-UTC.txt
  cat /sys/fs/cgroup/memory.peak > container-cumulative-memory-peak-bytes.txt
  echo "$name mesh generation and inspections complete"
done
date -u +%FT%TZ > /work/MESHES_COMPLETE.txt
