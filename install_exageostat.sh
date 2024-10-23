#!/bin/bash
build_nlopt=1
build_gsl=1
build_hwloc=1
build_starpu=1
build_chameleon=1
build_dplasma=1
# Run everything on Shaheen cdl5 login node - INTEL compiler

prgenv=${1:-cray-mkl}
case ${prgenv} in
intel)
  module swap PrgEnv-cray PrgEnv-intel
  #module load oneapi
  ;;
gnu|gnu-mkl)
  module swap PrgEnv-cray PrgEnv-gnu
  module swap gcc gcc/11.2.0
  module load mkl
  ;;
cray-mkl)
  module load gcc
  module unload cray-libsci
  #module load mkl
  export MKLROOT=/scratch/hpe/fthomas/css/.local/mkl/2019.5
  ;;
cray|cray-libsci)
  :
  ;;
aocc-aocl)
  module swap PrgEnv-cray PrgEnv-aocc
  module load gcc
  module unload aocc
  module load $CSS/modules/modulefiles/aocc/4.1.0
  module load $CSS/modules/modulefiles/aocl/4.1.0
  ;;
*)
  echo "Please specify a valid prgenv"
  exit 0
  ;;
esac
module list
export LC_ALL=en_US.UTF-8
export CRAYPE_LINK_TYPE=dynamic
export YACC="yacc -d"

#==================================================
mkdir codes
cd codes
export SETUP_DIR=$PWD
rm -rf exageostatr
cd $SETUP_DIR
> $SETUP_DIR/pkg_config.sh
case ${prgenv} in
intel)
  echo 'module swap PrgEnv-cray PrgEnv-intel' >>  $SETUP_DIR/pkg_config.sh
  echo 'module load openapi' >>  $SETUP_DIR/pkg_config.sh
  ;;
gnu|gnu-mkl)
  echo 'module swap PrgEnv-cray PrgEnv-gnu' >>  $SETUP_DIR/pkg_config.sh
  echo 'module swap gcc gcc/11.2.0' >>  $SETUP_DIR/pkg_config.sh
  echo 'module load mkl' >>  $SETUP_DIR/pkg_config.sh
  ;;
cray-mkl)
  echo "export MKLROOT=$MKLROOT" >>  $SETUP_DIR/pkg_config.sh
  ;;
aocc-aocl)
  echo 'module swap PrgEnv-cray PrgEnv-aocc' >>  $SETUP_DIR/pkg_config.sh
  echo 'module load aocl' >>  $SETUP_DIR/pkg_config.sh
  ;;
esac
echo 'module load cmake' >>  $SETUP_DIR/pkg_config.sh
echo 'export LC_ALL=en_US.UTF-8' >>  $SETUP_DIR/pkg_config.sh
echo 'export CRAYPE_LINK_TYPE=dynamic' >>  $SETUP_DIR/pkg_config.sh
#================================
banner nlopt
if [ $build_nlopt -eq 1 ] ; then
  if [ ! -f "nlopt-2.4.2.tar.gz" ]; then
    wget http://ab-initio.mit.edu/nlopt/nlopt-2.4.2.tar.gz
  fi
  rm -rf nlopt-2.4.2
  tar -zxvf nlopt-2.4.2.tar.gz
  cd nlopt-2.4.2
  [[ ! -d nlopt_install ]] || mkdir nlopt_install
  CC=gcc ./configure  --prefix=$PWD/nlopt_install --enable-shared --without-guile
  make -j 32
  make -j install
else
  cd nlopt-2.4.2
fi
export NLOPTROOT=$PWD
export PKG_CONFIG_PATH=$NLOPTROOT/nlopt_install/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$NLOPTROOT/nlopt_install/lib:$LD_LIBRARY_PATH
echo 'export PKG_CONFIG_PATH='$NLOPTROOT'/nlopt_install/lib/pkgconfig:$PKG_CONFIG_PATH' >> $SETUP_DIR/pkg_config.sh
echo 'export LD_LIBRARY_PATH='$NLOPTROOT'/nlopt_install/lib:$LD_LIBRARY_PATH' >> $SETUP_DIR/pkg_config.sh
#================================
banner gsl
cd $SETUP_DIR
if [ $build_gsl -eq 1 ] ; then
  if [ ! -f "gsl-2.4.tar.gz" ]; then
    wget https://ftp.gnu.org/gnu/gsl/gsl-2.4.tar.gz
  fi
  rm -rf gsl-2.4
  tar -zxvf gsl-2.4.tar.gz
  cd gsl-2.4
  [[ ! -d gsl_install ]] || mkdir gsl_install
  ./configure  --prefix=$PWD/gsl_install/
  make -j 96
  make -j install
else
  cd gsl-2.4
fi
GSLROOT=$PWD
export PKG_CONFIG_PATH=$GSLROOT/gsl_install/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$GSLROOT/gsl_install/lib:$LD_LIBRARY_PATH
echo 'export PKG_CONFIG_PATH='$GSLROOT'/gsl_install/lib/pkgconfig:$PKG_CONFIG_PATH' >> $SETUP_DIR/pkg_config.sh
echo 'export LD_LIBRARY_PATH='$GSLROOT'/gsl_install/lib:$LD_LIBRARY_PATH' >> $SETUP_DIR/pkg_config.sh
#================================
banner hwloc
cd $SETUP_DIR
if [ $build_hwloc -eq 1 ] ; then
  if [  ! -f "hwloc-2.5.0.tar.gz" ]; then
    wget https://download.open-mpi.org/release/hwloc/v2.5/hwloc-2.5.0.tar.gz
  fi
  rm -rf hwloc-2.5.0
  tar -zxvf  hwloc-2.5.0.tar.gz
  cd hwloc-2.5.0
  [[ ! -d hwloc_install ]] || mkdir hwloc_install
  ./configure  --prefix=$PWD/hwloc_install --disable-libxml2 -disable-pci --enable-shared=yes
  make -j 32
  make -j install
else
  cd hwloc-2.5.0
fi
export HWLOCROOT=$PWD
export PKG_CONFIG_PATH=$HWLOCROOT/hwloc_install/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$HWLOCROOT/hwloc_install/lib:$LD_LIBRARY_PATH
echo 'export PKG_CONFIG_PATH='$HWLOCROOT'/hwloc_install/lib/pkgconfig:$PKG_CONFIG_PATH' >> $SETUP_DIR/pkg_config.sh
echo 'export LD_LIBRARY_PATH='$HWLOCROOT'/hwloc_install/lib:$LD_LIBRARY_PATH' >> $SETUP_DIR/pkg_config.sh
#================================
banner starpu
cd $SETUP_DIR
if [ $build_starpu -eq 1 ] ; then
  if [ ! -f "starpu-1.2.9.tar.gz" ]; then
    wget https://files.inria.fr/starpu/starpu-1.2.9/starpu-1.2.9.tar.gz
  fi
  rm -rf starpu-1.2.9
  tar -zxvf starpu-1.2.9.tar.gz
  cd starpu-1.2.9
  perl -i -p -e "s/, struct starpu_sched_component \* to STARPU_ATTRIBUTE_UNUSED//" src/sched_policies/component_perfmodel_select.c
  perl -p -e "s/                function fstarpu_worker_get_by_devid/                integer function fstarpu_worker_get_by_devid/" -i include/fstarpu_mod.f90
  [[ ! -d starpu_install ]] || mkdir starpu_install
  #CFLAGS=-fPIC CXXFLAGS=-fPIC 
  perl -p -e "s/MPICC_LDFLAGS=/#MPICC_LDFLAGS=/" -i configure
  perl -p -e "s/if text x/if test x/" -i configure
  ./configure  --prefix=$PWD/starpu_install --disable-cuda --disable-opencl --enable-shared --disable-build-doc --disable-export-dynamic --disable-mpi-check  --with-mpicc=$(which cc) --with-mpifort=$(which ftn) --with-mpiexec=$(which srun) --disable-starpu-top --disable-full-gdb-information --enable-maxcpus=128 --enable-maxnodes=8
  module load cray-fftw
  case $prgenv in
  *aocc*) make -j 32 CC=cc LDFLAGS="-Wl,-allow-multiple-definition -L/opt/cray/pe/gcc/12.2.0/snos/lib64";;
  *) make -j 32 LDFLAGS="-Wl,-allow-multiple-definition";;
  esac
  make -j install
  module unload cray-fftw
else
  cd starpu-1.2.9
fi
export STARPUROOT=$PWD
export PKG_CONFIG_PATH=$STARPUROOT/starpu_install/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$STARPUROOT/starpu_install/lib:$LD_LIBRARY_PATH
export CPATH=$STARPUROOT/starpu_install/include:$CPATH
export CPATH=$STARPUROOT/starpu_install/include/starpu/1.2:$CPATH
echo 'export PKG_CONFIG_PATH='$STARPUROOT'/starpu_install/lib/pkgconfig:$PKG_CONFIG_PATH' >> $SETUP_DIR/pkg_config.sh
echo 'export LD_LIBRARY_PATH='$STARPUROOT'/starpu_install/lib:$LD_LIBRARY_PATH' >> $SETUP_DIR/pkg_config.sh
echo 'export CPATH='$STARPUROOT'/starpu_install/include:$CPATH' >> $SETUP_DIR/pkg_config.sh
echo 'export CPATH='$STARPUROOT'/starpu_install/include/starpu/1.2:$CPATH' >> $SETUP_DIR/pkg_config.sh
#************************************************************************ Install Chameleon - Stars-H - HiCMA
#module load intel
cd $SETUP_DIR
cd ../exageostat_parsec
#git submodule update --init --recursive

export EXAGEOSTATDEVDIR=$PWD
cd $EXAGEOSTATDEVDIR

export HICMADIR=$EXAGEOSTATDEVDIR/hicma
export CHAMELEONDIR=$EXAGEOSTATDEVDIR/hicma/chameleon
export STARSHDIR=$EXAGEOSTATDEVDIR/stars-h

case $prgenv in
intel)
  SCILIBS="-Wl,--no-as-needed;-L${MKLROOT}/lib;-lmkl_intel_lp64;-lmkl_core;-lmkl_intel_thread;-lpthread;-lm;-ldl"
  SCILIBS="-Wl,--no-as-needed;-L${MKLROOT}/lib;-lmkl_intel_lp64;-lmkl_core;-lmkl_sequential;-lpthread;-lm;-ldl"
  ;;
gnu|gnu-mkl|cray-mkl)
  SCILIBS="-Wl,--no-as-needed;-L${MKLROOT}/lib;-lmkl_gf_lp64;-lmkl_core;-lmkl_gnu_thread;-lpthread;-lm;-ldl"
  SCILIBS="-Wl,--no-as-needed;-L${MKLROOT}/lib;-lmkl_gf_lp64;-lmkl_core;-lmkl_sequential;-lpthread;-lm;-ldl"
  ;;
cray|cray-libsci)
  SCILIBS="-Wl,--no-as-needed;-L${CRAY_LIBSCI_DIR}/cray/12.0/x86_64/lib;-lsci_cray_mp;-lpthread;-lm;-ldl"
  ;;
aocc-aocl)
  SCILIBS="-Wl,--no-as-needed;-L${AOCL_ROOT}/lib;-lblis;-lflame;-ltmg;-lmissing;-lpthread;-lm;-ldl"
  ;;
esac

banner CHAMELEON
## CHAMELEON
cd $CHAMELEONDIR
if [ $build_chameleon -eq 1 ] ; then
  rm -rf build
  mkdir -p build/install_dir
  cd build

  case $prgenv in
  intel)
  LDFLAGS="-Wl,-allow-shlib-undefined -lrt" cmake .. -DCMAKE_Fortran_FLAGS="-nofor-main"  -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_FLAGS='-g' -DCMAKE_C_FLAGS='-g' -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_INSTALL_PREFIX=$PWD/install_dir -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=ON -DCHAMELEON_ENABLE_EXAMPLE=ON -DCHAMELEON_ENABLE_TESTING=ON -DCHAMELEON_ENABLE_TIMING=ON -DCHAMELEON_USE_MPI=ON -DCHAMELEON_USE_CUDA=OFF -DCHAMELEON_USE_MAGMA=OFF -DCHAMELEON_SCHED_QUARK=OFF -DCHAMELEON_SCHED_STARPU=ON -DCHAMELEON_USE_FXT=OFF -DSTARPU_DIR=$STARPUROOT/starpu_install -DBLAS_LIBRARIES="$SCILIB" -DBLAS_COMPILER_FLAGS="-m64;-I${MKLROOT}/include" -DLAPACK_LIBRARIES="$SCILIBS" -DCBLAS_LIBRARIES="$SCILIBS" -DLAPACKE_DIR="${MKLROOT}" -DTMG_DIR="${MKLROOT}" -DMORSE_VERBOSE_FIND_PACKAGE=ON -DMPI_C_COMPILER=$(which cc) -DBUILD_SHARED_LIBS=ON -DPYTHON_EXECUTABLE=/usr/bin/python3.6
  ;;
  gnu|gnu-mkl|cray-mkl)
  LDFLAGS="-Wl,-allow-shlib-undefined -lrt" cmake .. -DCMAKE_Fortran_FLAGS="-O2"  -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_FLAGS="-g -I${MKLROOT}/include" -DCMAKE_C_FLAGS="-g -I${MKLROOT}/include" -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_INSTALL_PREFIX=$PWD/install_dir -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=ON -DCHAMELEON_ENABLE_EXAMPLE=ON -DCHAMELEON_ENABLE_TESTING=ON -DCHAMELEON_ENABLE_TIMING=ON -DCHAMELEON_USE_MPI=ON -DCHAMELEON_USE_CUDA=OFF -DCHAMELEON_USE_MAGMA=OFF -DCHAMELEON_SCHED_QUARK=OFF -DCHAMELEON_SCHED_STARPU=ON -DCHAMELEON_USE_FXT=OFF -DSTARPU_DIR=$STARPUROOT/starpu_install -DBLAS_LIBRARIES="$SCILIBS" -DBLAS_COMPILER_FLAGS="-m64;-I${MKLROOT}/include" -DLAPACK_LIBRARIES="$SCILIBS" -DCBLAS_LIBRARIES="$SCILIBS" -DLAPACKE_DIR="${MKLROOT}" -DTMG_DIR="${MKLROOT}" -DMORSE_VERBOSE_FIND_PACKAGE=ON -DMPI_C_COMPILER=$(which cc) -DBUILD_SHARED_LIBS=ON -DPYTHON_EXECUTABLE=/usr/bin/python3.6
  ;;
  aocc-aocl)
  LDFLAGS="-L/opt/cray/pe/gcc/12.2.0/snos/lib64 -Wl,-allow-shlib-undefined -lrt" cmake .. -DCMAKE_Fortran_FLAGS="-O2"  -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_FLAGS="-g -I${AOCL_ROOT}/include" -DCMAKE_C_FLAGS="-g -I${AOCL_ROOT}/include" -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_INSTALL_PREFIX=$PWD/install_dir -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=ON -DCHAMELEON_ENABLE_EXAMPLE=ON -DCHAMELEON_ENABLE_TESTING=ON -DCHAMELEON_ENABLE_TIMING=ON -DCHAMELEON_USE_MPI=ON -DCHAMELEON_USE_CUDA=OFF -DCHAMELEON_USE_MAGMA=OFF -DCHAMELEON_SCHED_QUARK=OFF -DCHAMELEON_SCHED_STARPU=ON -DCHAMELEON_USE_FXT=OFF -DSTARPU_DIR=$STARPUROOT/starpu_install -DBLAS_LIBRARIES="$SCILIBS" -DBLAS_COMPILER_FLAGS="-m64;-I${AOCL_ROOT}/include" -DLAPACK_LIBRARIES="$SCILIBS" -DCBLAS_LIBRARIES="$SCILIBS" -DLAPACKE_DIR="${AOCL_ROOT}" -DTMG_DIR="${AOCL_ROOT}" -DMORSE_VERBOSE_FIND_PACKAGE=ON -DMPI_C_COMPILER=$(which cc) -DBUILD_SHARED_LIBS=ON -DPYTHON_EXECUTABLE=/usr/bin/python3.6
  ;;
  cray|cray-libsci)
  LDFLAGS=-lrt cmake .. -DCMAKE_Fortran_FLAGS="-O2"  -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_FLAGS='-g' -DCMAKE_C_FLAGS='-g' -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_INSTALL_PREFIX=$PWD/install_dir -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=ON -DCHAMELEON_ENABLE_EXAMPLE=ON -DCHAMELEON_ENABLE_TESTING=ON -DCHAMELEON_ENABLE_TIMING=ON -DCHAMELEON_USE_MPI=ON -DCHAMELEON_USE_CUDA=OFF -DCHAMELEON_USE_MAGMA=OFF -DCHAMELEON_SCHED_QUARK=OFF -DCHAMELEON_SCHED_STARPU=ON -DCHAMELEON_USE_FXT=OFF -DSTARPU_DIR=$STARPUROOT/starpu_install -DBLAS_LIBRARIES="$SCILIBS" -DBLAS_COMPILER_FLAGS="-m64;-I${CRAY_LIBSCI_DIR}/include" -DLAPACK_LIBRARIES="$SCILIBS" -DCBLAS_DIR="${CRAY_LIBSCI_DIR}/cray/12.0/x86_64" -DLAPACKE_DIR="${CRAY_LIBSCI_DIR}/cray/12.0/x86_64" -DTMG_DIR="${CRAY_LIBSCI_DIR}/cray/12.0/x86_64" -DMORSE_VERBOSE_FIND_PACKAGE=ON -DMPI_C_COMPILER=$(which cc) -DBUILD_SHARED_LIBS=ON -DPYTHON_EXECUTABLE=/usr/bin/python3.6
  ;;
  *)
  LDFLAGS=-lrt cmake .. -DCMAKE_Fortran_FLAGS="-O2"  -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_FLAGS="-g -I${MKLROOT}/include" -DCMAKE_C_FLAGS="-g -I${MKLROOT}/include" -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_INSTALL_PREFIX=$PWD/install_dir -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=ON -DCHAMELEON_ENABLE_EXAMPLE=ON -DCHAMELEON_ENABLE_TESTING=ON -DCHAMELEON_ENABLE_TIMING=ON -DCHAMELEON_USE_MPI=ON -DCHAMELEON_USE_CUDA=OFF -DCHAMELEON_USE_MAGMA=OFF -DCHAMELEON_SCHED_QUARK=OFF -DCHAMELEON_SCHED_STARPU=ON -DCHAMELEON_USE_FXT=OFF -DSTARPU_DIR=$STARPUROOT/starpu_install -DBLAS_LIBRARIES="-Wl,--no-as-needed;-L${MKLROOT}/lib;-lmkl_intel_lp64;-lmkl_core;-lmkl_gnu_thread;-lpthread;-lm;-ldl" -DBLAS_COMPILER_FLAGS="-m64;-I${MKLROOT}/include" -DLAPACK_LIBRARIES="-Wl,--no-as-needed;-L${MKLROOT}/lib;-lmkl_intel_lp64;-lmkl_core;-lmkl_gnu_thread;-lpthread;-lm;-ldl" -DCBLAS_DIR="${MKLROOT}" -DLAPACKE_DIR="${MKLROOT}" -DTMG_DIR="${MKLROOT}" -DMORSE_VERBOSE_FIND_PACKAGE=ON -DMPI_C_COMPILER=$(which cc) -DBUILD_SHARED_LIBS=ON -DPYTHON_EXECUTABLE=/usr/bin/python3.6
  ;;
  esac

  make -j 32 install # CHAMELEON parallel build seems to be fixed
  make install
fi

export PKG_CONFIG_PATH=$CHAMELEONDIR/build/install_dir/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=$CHAMELEONDIR/build/install_dir/lib/:$LD_LIBRARY_PATH
export CPATH=$CHAMELEONDIR/build/install_dir/include/coreblas:$CPATH
export PATH=$CHAMELEONDIR/build/install_dir/include/coreblas:$PATH

echo 'export PKG_CONFIG_PATH='$CHAMELEONDIR'/build/install_dir/lib/pkgconfig:$PKG_CONFIG_PATH' >> $SETUP_DIR/pkg_config.sh
echo 'export LD_LIBRARY_PATH='$CHAMELEONDIR'/build/install_dir/lib/:$LD_LIBRARY_PATH' >> $SETUP_DIR/pkg_config.sh
echo 'export CPATH='$CHAMELEONDIR'/build/install_dir/include/coreblas:$CPATH' >> $SETUP_DIR/pkg_config.sh
echo 'export PATH='$CHAMELEONDIR'/build/install_dir/include/coreblas:$PATH' >> $SETUP_DIR/pkg_config.sh

#==============================================
#module unload cmake
# cmake 3.18 
#export PATH=/lustre/project/k1205/lei/software/CMake/bin:$PATH

banner dplasma
cd $SETUP_DIR
cd ../exageostat_parsec/dplasma
if [ $build_dplasma -eq 1 ] ; then
  #git checkout 2fbef5d868f51f08300231045f1310e469821bb0
  #git submodule update --init --recursive
  rm -rf build
  mkdir build
  cd build
  case $prgenv in
  intel)
  #cmake .. -DCMAKE_Fortran_FLAGS="-nofor-main" -DDPLASMA_PRECISIONS="d" -DBUILD_SHARED_LIBS=ON -DPARSEC_WITH_HEADER_FILES=ON -DPARSEC_WITH_DEVEL_HEADERS=ON -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn  -DPARSEC_GPU_WITH_CUDA=OFF -DPARSEC_HAVE_CUDA=off -DCMAKE_BUILD_TYPE=Debug -DPARSEC_DEBUG_PARANOID=ON -DPARSEC_DEBUG_NOISIER=ON -DCMAKE_INSTALL_PREFIX=$PWD/install_dir
  LDFLAGS="-Wl,-allow-multiple-definition" cmake .. -DCMAKE_Fortran_FLAGS="-nofor-main" -DDPLASMA_PRECISIONS="d" -DBUILD_SHARED_LIBS=ON -DPARSEC_WITH_HEADER_FILES=ON -DPARSEC_WITH_DEVEL_HEADERS=ON -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn  -DPARSEC_GPU_WITH_CUDA=OFF -DPARSEC_HAVE_CUDA=off -DCMAKE_INSTALL_PREFIX=$PWD/install_dir -DPYTHON_EXECUTABLE=/usr/bin/python3.6
  ;;
  gnu|gnu-mkl|cray-mkl)
  LDFLAGS="-Wl,-allow-multiple-definition -Wl,-allow-shlib-undefined" cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_FLAGS="-O2" -DDPLASMA_PRECISIONS="d" -DBUILD_SHARED_LIBS=ON -DPARSEC_WITH_HEADER_FILES=ON -DPARSEC_WITH_DEVEL_HEADERS=ON -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn -DLAPACK_LIBRARIES="$SCILIBS" -DBLAS_LIBRARIES="$SCILIBS" -DCBLAS_LIBRARIES="$SCILIBS" -DLAPACKE_DIR="${MKLROOT}" -DLAPACKE_CBLAS_LIB="$SCILIBS" -DLAPACKE_LAPACKE_LIB="$SCILIBS" -DPARSEC_GPU_WITH_CUDA=OFF -DPARSEC_HAVE_CUDA=off -DPARSEC_PROF_GRAPHER=OFF -DCMAKE_INSTALL_PREFIX=$PWD/install_dir -DPYTHON_EXECUTABLE=/usr/bin/python3.6
  ;;
  aocc-aocl)
  LDFLAGS="-L/opt/cray/pe/gcc/12.2.0/snos/lib64 -Wl,-allow-multiple-definition -Wl,-allow-shlib-undefined" cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_FLAGS="-O2" -DDPLASMA_PRECISIONS="d" -DBUILD_SHARED_LIBS=ON -DPARSEC_WITH_HEADER_FILES=ON -DPARSEC_WITH_DEVEL_HEADERS=ON -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn -DLAPACK_LIBRARIES="$SCILIBS" -DBLAS_LIBRARIES="$SCILIBS" -DCBLAS_LIBRARIES="$SCILIBS" -DLAPACKE_DIR="${AOCL_ROOT}" -DLAPACKE_CBLAS_LIB="$SCILIBS" -DLAPACKE_LAPACKE_LIB="$SCILIBS" -DPARSEC_GPU_WITH_CUDA=OFF -DPARSEC_HAVE_CUDA=off -DPARSEC_PROF_GRAPHER=OFF -DCMAKE_INSTALL_PREFIX=$PWD/install_dir -DPYTHON_EXECUTABLE=/usr/bin/python3.6
  ;;
  *)
  LDFLAGS="-Wl,-allow-multiple-definition" cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_FLAGS="-O2" -DDPLASMA_PRECISIONS="d" -DBUILD_SHARED_LIBS=ON -DPARSEC_WITH_HEADER_FILES=ON -DPARSEC_WITH_DEVEL_HEADERS=ON -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn -DLAPACK_LIBRARIES="$SCILIBS" -DCBLAS_DIR="${MKLROOT}" -DLAPACKE_DIR="${MKLROOT}" -DLAPACKE_CBLAS_LIB="$SCILIBS" -DLAPACKE_LAPACKE_LIB="$SCILIBS" -DPARSEC_GPU_WITH_CUDA=OFF -DPARSEC_HAVE_CUDA=off -DCMAKE_INSTALL_PREFIX=$PWD/install_dir -DPYTHON_EXECUTABLE=/usr/bin/python3.6
  ;;
  esac

  make -j 32 VERBOSE=1
  make install
  cd ..
fi

export PKG_CONFIG_PATH=$PWD/build/install_dir/lib/pkgconfig:$PWD/build/install_dir/lib64/pkgconfig:$PKG_CONFIG_PATH
export PaRSEC_ROOT=$PWD/build/install_dir
echo "export PKG_CONFIG_PATH=$PWD/build/install_dir/lib/pkgconfig:\$PKG_CONFIG_PATH"  >> $SETUP_DIR/pkg_config.sh
echo "export PKG_CONFIG_PATH=$PWD/build/install_dir/lib64/pkgconfig:\$PKG_CONFIG_PATH"  >> $SETUP_DIR/pkg_config.sh
echo "export PaRSEC_ROOT=$PWD/build/install_dir"  >> $SETUP_DIR/pkg_config.sh
#===================================
banner ExaGeoStat
cd $SETUP_DIR
#ExaGeoStat
export PATH=$CHAMELEONDIR/build/install_dir/include/coreblas:$PATH
export CPATH=$CHAMELEONDIR/build/install_dir/include/coreblas:$CPATH

echo 'export PATH='$CHAMELEONDIR'/build/install_dir/include/coreblas:$PATH'  >> $SETUP_DIR/pkg_config.sh
echo 'export CPATH='$CHAMELEONDIR'/build/install_dir/include/coreblas:$CPATH'  >> $SETUP_DIR/pkg_config.sh
cd $EXAGEOSTATDEVDIR
rm -rf build
mkdir build
cd build
case $prgenv in
intel)
cmake ..  -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_FLAGS='-g' -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_INSTALL_PREFIX=$PWD/install_dir -DEXAGEOSTAT_USE_MPI=OFF -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=ON    -DEXAGEOSTAT_USE_HICMA=OFF -DEXAGEOSTAT_USE_NETCDF=OFF -DLAPACKE_INCLUDE_DIRS="-m64;-I${MKLROOT}/include"  -DBLAS_LIBRARIES="$SCILIBS" -DBLAS_COMPILER_FLAGS="-m64;-I${MKLROOT}/include" -DLAPACK_LIBRARIES="$SCILIBS" -DCBLAS_DIR="${MKLROOT}" -DLAPACKE_DIR="${MKLROOT}" -DTMG_DIR="${MKLROOT}" -DEXAGEOSTAT_USE_DPLASMA=ON -DLAPACKE_LIBRARIES="$SCILIBS"
;;
cray|cray-libsci)
LDFLAGS="-Wl,-allow-multiple-definition" cmake ..  -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_FLAGS='-g' -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_INSTALL_PREFIX=$PWD/install_dir -DEXAGEOSTAT_USE_MPI=OFF -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=ON    -DEXAGEOSTAT_USE_HICMA=OFF -DEXAGEOSTAT_USE_NETCDF=OFF -DLAPACKE_INCLUDE_DIRS="-m64;-I${CRAY_LIBSCI_DIR}/cray/12.0/x86_64/include"  -DBLAS_LIBRARIES="$SCILIBS" -DBLAS_COMPILER_FLAGS="-m64;-I${CRAY_LIBSCI_DIR}/cray/12.0/x86_64/include" -DLAPACK_LIBRARIES="$SCILIBS" -DCBLAS_DIR="${CRAY_LIBSCI_DIR}/cray/12.0/x86_64" -DLAPACKE_DIR="${CRAY_LIBSCI_DIR}/cray/12.0/x86_64" -DTMG_DIR="${CRAY_LIBSCI_DIR}/cray/12.0/x86_64" -DEXAGEOSTAT_USE_DPLASMA=ON -DLAPACKE_LIBRARIES="$SCILIBS"
for i in $(grep -rlw "lstarpu-1.2" . --include=link.txt);do
  perl -p -e "s/lstarpu-1.2/lstarpumpi-1.2 -lstarpu-1.2/" -i $i
done
;;
gnu|gnu-mkl|cray-mkl)
LDFLAGS="-Wl,-allow-multiple-definition -Wl,-allow-shlib-undefined" cmake ..  -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_CPP_FLAGS="-DKJT_PROFILER" -DCMAKE_CXX_FLAGS='-g' -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_INSTALL_PREFIX=$PWD/install_dir -DEXAGEOSTAT_USE_MPI=OFF -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=ON    -DEXAGEOSTAT_USE_HICMA=OFF -DEXAGEOSTAT_USE_NETCDF=OFF -DLAPACKE_INCLUDE_DIRS="-m64;-I${MKLROOT}/include"  -DBLAS_LIBRARIES="$SCILIBS" -DBLAS_COMPILER_FLAGS="-m64;-I${MKLROOT}/include" -DLAPACK_LIBRARIES="$SCILIBS" -DCBLAS_LIBRARIES="$SCILIBS" -DLAPACKE_DIR="${MKLROOT}" -DTMG_DIR="${MKLROOT}" -DEXAGEOSTAT_USE_DPLASMA=ON -DLAPACKE_LIBRARIES="$SCILIBS" -DPYTHON_EXECUTABLE=/usr/bin/python3.6 -DLAPACKE_CBLAS_LIB="$SCILIBS" -DLAPACKE_LAPACKE_LIB="$SCILIBS" -DHWLOC_DIR="$HWLOCROOT/hwloc_install"
for i in $(grep -rlw "lstarpu-1.2" . --include=link.txt);do
  perl -p -e "s/lstarpu-1.2/lstarpumpi-1.2 -lstarpu-1.2/" -i $i
done
;;
aocc-aocl)
LDFLAGS="-L/opt/cray/pe/gcc/12.2.0/snos/lib64 -Wl,-allow-multiple-definition -Wl,-allow-shlib-undefined" cmake ..  -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_CPP_FLAGS="-DKJT_PROFILER" -DCMAKE_CXX_FLAGS='-g' -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_INSTALL_PREFIX=$PWD/install_dir -DEXAGEOSTAT_USE_MPI=OFF -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=ON    -DEXAGEOSTAT_USE_HICMA=OFF -DEXAGEOSTAT_USE_NETCDF=OFF -DLAPACKE_INCLUDE_DIRS="-m64;-I${AOCL_ROOT}/include"  -DBLAS_LIBRARIES="$SCILIBS" -DBLAS_COMPILER_FLAGS="-m64;-I${AOCL_ROOT}/include" -DLAPACK_LIBRARIES="$SCILIBS" -DCBLAS_LIBRARIES="$SCILIBS" -DLAPACKE_DIR="${AOCL_ROOT}" -DTMG_DIR="${AOCL_ROOT}" -DEXAGEOSTAT_USE_DPLASMA=ON -DLAPACKE_LIBRARIES="$SCILIBS" -DPYTHON_EXECUTABLE=/usr/bin/python3.6 -DLAPACKE_CBLAS_LIB="$SCILIBS" -DLAPACKE_LAPACKE_LIB="$SCILIBS" -DHWLOC_DIR="$HWLOCROOT/hwloc_install"
for i in $(grep -rlw "lstarpu-1.2" . --include=link.txt);do
  perl -p -e "s/lstarpu-1.2/lstarpumpi-1.2 -lstarpu-1.2/" -i $i
done
;;
*)
LDFLAGS="-Wl,-allow-multiple-definition" cmake ..  -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_FLAGS='-g' -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_INSTALL_PREFIX=$PWD/install_dir -DEXAGEOSTAT_USE_MPI=OFF -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=ON    -DEXAGEOSTAT_USE_HICMA=OFF -DEXAGEOSTAT_USE_NETCDF=OFF -DLAPACKE_INCLUDE_DIRS="-m64;-I${MKLROOT}/include"  -DBLAS_LIBRARIES="$SCILIBS" -DBLAS_COMPILER_FLAGS="-m64;-I${MKLROOT}/include" -DLAPACK_LIBRARIES="$SCILIBS" -DCBLAS_DIR="${MKLROOT}" -DLAPACKE_DIR="${MKLROOT}" -DTMG_DIR="${MKLROOT}" -DEXAGEOSTAT_USE_DPLASMA=ON -DLAPACKE_LIBRARIES="$SCILIBS"
for i in $(grep -rlw "lstarpu-1.2" . --include=link.txt);do
  perl -p -e "s/lstarpu-1.2/lstarpumpi-1.2 -lstarpu-1.2/" -i $i
done
;;
esac
make -j 32 VERBOSE=1
make install
