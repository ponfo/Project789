#!/bin/bash

rm $2/projectData.dat

cat <<EOT>> projectData.dat
$1
$3
EOT


./../../main $2/$1
