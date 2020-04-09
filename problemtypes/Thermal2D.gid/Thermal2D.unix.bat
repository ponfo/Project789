#!/bin/bash

rm $2/projectData.dat

cat <<EOT>> projectData.dat
$1
$3
EOT

$3/../../applications/Thermal2DGid/main $2/$1
