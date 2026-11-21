#!/bin/bash
set --
source /usr/lib/openfoam/openfoam2412/etc/bashrc
set -euo pipefail
for name in B00 A01; do
  cd /work/"$name"
  test ! -e log.snappyHexMesh
  surfaceCheck constant/triSurface/volume.stl > log.surfaceCheck 2>&1
  grep -q 'Surface is closed' log.surfaceCheck
  surfaceFeatureExtract > log.surfaceFeatureExtract 2>&1
  blockMesh > log.blockMesh 2>&1
  decomposePar -no-fields > log.decomposePar-mesh 2>&1
  mpirun --mca btl_vader_single_copy_mechanism none -np 4 snappyHexMesh -parallel -overwrite > log.snappyHexMesh 2>&1
  reconstructParMesh -constant -noZero > log.reconstructParMesh 2>&1
  checkMesh -allTopology > log.checkMesh 2>&1
  checkMesh -allTopology -allGeometry > log.checkMesh-additional 2>&1 || true
  grep -q 'Mesh OK' log.checkMesh
  echo "$name mesh complete"
done
