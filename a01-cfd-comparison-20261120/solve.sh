#!/bin/bash
set --
source /usr/lib/openfoam/openfoam2412/etc/bashrc
set -euo pipefail
for name in B00 A01; do
  cd /work/"$name"
  grep -q 'Mesh OK' log.checkMesh
  test ! -e log.simpleFoam
  cp system/fvSchemes system/fvSchemes-second-order
  cp system/controlDict system/controlDict-production
  cp system/fvSolution system/fvSolution-production
  decomposePar -force > log.decomposePar-solver 2>&1
  mpirun --mca btl_vader_single_copy_mechanism none -np 4 potentialFoam -parallel -initialiseUBCs > log.potentialFoam 2>&1
  # Matched first-order startup for both geometries, then common second order.
  sed -i 's/bounded Gauss linearUpwind grad(U)/bounded Gauss upwind/' system/fvSchemes
  sed -i 's/residualControl {p 1e-4; U 1e-5; k 1e-5; omega 1e-5;}//' system/fvSolution
  foamDictionary system/controlDict -entry endTime -set 100 > log.control-warmup 2>&1
  mpirun --mca btl_vader_single_copy_mechanism none -np 4 simpleFoam -parallel > log.simpleFoam-warmup 2>&1
  cp system/fvSchemes-second-order system/fvSchemes
  cp system/controlDict-production system/controlDict
  cp system/fvSolution-production system/fvSolution
  foamDictionary system/controlDict -entry startFrom -set latestTime > log.control-production 2>&1
  mpirun --mca btl_vader_single_copy_mechanism none -np 4 simpleFoam -parallel > log.simpleFoam 2>&1
  reconstructPar -latestTime > log.reconstructPar 2>&1
  foamToVTK -latestTime -fields '(U p)' > log.foamToVTK 2>&1
  previous=$(find processor0 -maxdepth 1 -type d -printf '%f\n' | awk '/^[0-9]+$/ && $0>0' | sort -n | tail -2 | head -1)
  latest=$(find processor0 -maxdepth 1 -type d -printf '%f\n' | awk '/^[0-9]+$/ && $0>0' | sort -n | tail -1)
  if test "$previous" != "$latest"; then
    reconstructPar -time "$previous" -fields '(U p)' > log.reconstructPrevious 2>&1
    foamToVTK -time "$previous" -fields '(U p)' > log.foamToVTK-previous 2>&1
  fi
  echo "$name solve complete"
done
