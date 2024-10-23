#!/bin/bash

# custom environment used during CI tests with gitlab ci

# these paths may depend on the runner used, please be careful and add
# the necessary if blocks depending on the machine

. /etc/profile.d/spack.sh

export QUARK_DIR=/builds/install/quark
export PARSEC_DIR=/builds/install/parsec
export STARPU_DIR=/builds/install/starpu
export STARPU_SILENT=1

if [ "$1" == "simu" ]; then
  export STARPU_DIR=/builds/install/starpu-simgrid
  export SIMGRID_DIR=/builds/install/simgrid
fi

export PKG_CONFIG_PATH=$PARSEC_DIR/lib/pkgconfig:$PKG_CONFIG_PATH
export PKG_CONFIG_PATH=$STARPU_DIR/lib/pkgconfig:$PKG_CONFIG_PATH
