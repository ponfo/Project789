#!/bin/bash

rm $2/EULER.DAT

cat <<EOT>> EULER.DAT
$1
EOT

$3/../../main $2/$1
