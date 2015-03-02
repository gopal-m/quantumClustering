#!/bin/bash

# This script requires the devtools package in order to function!

VERSION=$(grep 'Version' DESCRIPTION | cut -d' ' -f2)
FILENAME=quantumClustering_${VERSION}.tar.gz
DIR=$(pwd)

echo "library(devtools); devtools::document('${DIR}'); devtools::build('${DIR}'); install.packages('${DIR}/../${FILENAME}')" | R --no-save --slave
