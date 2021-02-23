
. ./adaptation_parameters.sh

. ./spack-sgpp/share/spack/setup-env.sh
spack load sgpp@adaptive_combigrid_convenience

python3 initialize_adaptive_combigrid.py