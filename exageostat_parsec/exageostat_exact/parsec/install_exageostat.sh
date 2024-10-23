#!/bin/bash

DEBUG=1
if [ $DEBUG -eq 1 ]; then
    BUILD_TYPE=Debug
else
    BUILD_TYPE=Release
fi


STARPU=starpu-1.2.5
EXAGEOSTAT=exageostat_parsec
CHAMELEON=chameleon
HICMA=hicma
STARS=stars-h
DPLASMA=dplasma

HOME=/home/qcao3
INSTALL=$HOME/opt
home=$HOME/lab/$EXAGEOSTAT

# install spack
git clone https://github.com/spack/spack
echo source $HOME/spack/share/spack/setup-env.sh | tee -a $HOME/.bashrc
source $HOME/.bashrc

# istall cmake
spack install cmake 
echo "spack load cmake" | tee -a $HOME/.bashrc
spack load cmake 

# install openblas
#spack install openblas 
#echo "spack load openblas" | tee -a $HOME/.bashrc
#spack load openblas 

# install hwloc 
#spack install hwloc 
#echo "spack load hwloc" | tee -a $HOME/.bashrc
#spack load hwloc 

#install nlopt
spack install nlopt 
echo "spack load nlopt" | tee -a $HOME/.bashrc
spack load nlopt

#install gsl
spack install gsl
echo "spack load gsl" | tee -a $HOME/.bashrc
spack load gsl 

#install netcdf-c 
spack install netcdf-c 
echo "spack load netcdf-c" | tee -a $HOME/.bashrc
spack load netcdf-c 

# Install starpu
cd $home/../
mkdir -p $INSTALL/$STARPU
wget http://starpu.gforge.inria.fr/files/$STARPU/$STARPU".tar.gz"
tar -zxvf $STARPU".tar.gz"
cd starpu-1.2.5
./configure --disable-cuda --disable-opencl --prefix=$INSTALL/$STARPU
make -j 8 && make install
export PKG_CONFIG_PATH=$INSTALL/$STARPU/lib/pkgconfig:$PKG_CONFIG_PATH
echo "export PKG_CONFIG_PATH=$INSTALL/$STARPU/lib/pkgconfig:\$PKG_CONFIG_PATH" | tee -a $HOME/.bashrc

# updata submodule 
cd $home
git submodule update --init --recursive

# Install chameleon 
mkdir -p $INSTALL/$CHAMELEON
cd $home
cd $HICMA 
cd $CHAMELEON
mkdir -p build && cd build
cmake .. -DCHAMELEON_USE_MPI=OFF -DCHAMELEON_ENABLE_EXAMPLE=OFF -DCHAMELEON_ENABLE_TESTING=OFF -DCHAMELEON_ENABLE_TIMING=OFF -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=$INSTALL/$CHAMELEON
make -j 8 && make install
export PKG_CONFIG_PATH=$INSTALL/$CHAMELEON/lib/pkgconfig:$PKG_CONFIG_PATH
echo "export PKG_CONFIG_PATH=$INSTALL/$CHAMELEON/lib/pkgconfig:\$PKG_CONFIG_PATH" | tee -a $HOME/.bashrc

# Install star-h 
mkdir -p $INSTALL/$STARS
cd $home
cd $STARS
mkdir -p build && cd build
cmake .. -DCMAKE_C_FLAGS=-fPIC -DEXAMPLES=OFF -DTESTING=OFF -DMPI=OFF -DCMAKE_INSTALL_PREFIX=$INSTALL/$STARS
make -j 8 && make install
export PKG_CONFIG_PATH=$INSTALL/$STARS/lib/pkgconfig:$PKG_CONFIG_PATH
echo "export PKG_CONFIG_PATH=$INSTALL/$STARS/lib/pkgconfig:\$PKG_CONFIG_PATH" | tee -a $HOME/.bashrc

# Install hicma 
mkdir -p $INSTALL/$HICMA
cd $home
cd $HICMA
mkdir -p build && cd build
cmake .. -DBUILD_SHARED_LIBS=ON -DHICMA_ENABLE_TESTING=OFF -DHICMA_ENABLE_TIMING=OFF -DCMAKE_INSTALL_PREFIX=$INSTALL/$HICMA
make -j 8 && make install
export PKG_CONFIG_PATH=$INSTALL/$HICMA/lib/pkgconfig:$PKG_CONFIG_PATH
echo "export PKG_CONFIG_PATH=$INSTALL/$HICMA/lib/pkgconfig:\$PKG_CONFIG_PATH" | tee -a $HOME/.bashrc

# Install dplasma 
# commit for dplasma: 2fbef5d868f51f08300231045f1310e469821bb0
# commit for parsec:  f448893293f112d24e269e23137f384d40e973cc
mkdir -p $INSTALL/$DPLASMA
cd $home
cd $DPLASMA
git submodule update --init --recursive
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DDPLASMA_PRECISIONS="s;d" -DBLA_VENDOR=Intel10_64lp_seq -DCMAKE_INSTALL_PREFIX=$INSTALL/$DPLASMA
make -j 8 && make install
export PKG_CONFIG_PATH=$INSTALL/$DPLASMA/lib/pkgconfig:$PKG_CONFIG_PATH
echo "export PKG_CONFIG_PATH=$INSTALL/$DPLASMA/lib/pkgconfig:\$PKG_CONFIG_PATH" | tee -a $HOME/.bashrc
echo "export PKG_CONFIG_PATH=$INSTALL/$DPLASMA/dplasma/lib/pkgconfig:\$PKG_CONFIG_PATH" | tee -a $HOME/.bashrc
echo "export PaRSEC_ROOT=$INSTALL/$DPLASMA" | tee -a $HOME/.bashrc

# Install exageostat 
mkdir -p $INSTALL/$EXAGEOSTAT
cd $home
mkdir -p build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/home/qcao3/opt/exageostat -DEXAGEOSTAT_SCHED_STARPU=ON   -DEXAGEOSTAT_USE_NETCDF=ON -DEXAGEOSTAT_USE_HICMA=ON -DEXAGEOSTAT_USE_DPLASMA=ON -DCMAKE_BUILD_TYPE=$BUILD_TYPE
make -j 8 && make install
export PKG_CONFIG_PATH=$INSTALL/$EXAGEOSTAT/lib/pkgconfig:$PKG_CONFIG_PATH
echo "export PKG_CONFIG_PATH=$INSTALL/$EXAGEOSTAT/lib/pkgconfig:\$PKG_CONFIG_PATH" | tee -a $HOME/.bashrc

