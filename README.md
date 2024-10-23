# Build instructions

The [install_exageostat.sh](install_exageostat.sh), is our build script for ExaGeoStat_CPU. It is meant to be run from the top directory. It is currently set to build ExaGeoStat with the Cray compiler and Intel's MKL-2019.5. Make sure you adjust the MKLROOT variable near the top of the script. We have tried newer versions of MKL but they don't improve the performance. It builds nlopt, gsl, hwloc, starpu, chameleon and dplasma.

# Run instructions

The [launchall](launchall) job script will submit multiple ExaGeoStat_CPU jobs. You can tweak the tests bash array to select the jobs you wish to run. Please check MKLROOT near the top of the [job.template](job.template) file.

# Tips

  * The [job.template](job.template) script uses the [interleave.sh](interleave.sh) script to allow memory interleaving when an MPI rank uses cores belonging to multiple NUMA domains.
  * We also use a utility (makemask) to create CPU masks for the srun command. See ```/scratch/hpe/fthomas/makemask/``` on the system itself.
  * We have tried to set a few parsec parameters but it seems only the scheduler (PARSEC_SCHED) makes a difference.
  * We specify 5 parameters for every ExaGeoStat run:
    * the matrix size: 300K, 1.1M, 1.5M, 3.4M or 5.6M
    * the tile size: 320, 480 or whatever we found to work best
    * the number of nodes
   * the number of ranks per node: 1 being the best value for us
   * the depop parameter, which specifies how many cores are left unused. 1 is a minimum, for running the MPI thread, but we found that leaving some cores unused was beneficial for small matrices.
  * The job_* scripts sourced in the [job.template](job.template) script reside under ~x_thomasf
