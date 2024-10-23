#/bin/bash

shopt -s extglob

rm -rf codes/gsl-2.4 codes/hwloc-2.5.0 codes/nlopt-2.4.2 codes/starpu-1.2.9
rm -rf ./exageostat_parsec/hicma/chameleon/build
rm -rf ./exageostat_parsec/dplasma/build
rm -rf ./exageostat_parsec/build
rm -f 4.ExaGeoStat_CPU.+([0-9])
rm -rf shaheen3/runs/+([0-9])
rm -rf runs/+([0-9])

