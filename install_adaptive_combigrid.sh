#!/bin/bash

git clone https://github.com/spack/spack.git spack-sgpp
cd spack-sgpp

git apply ../adaptive_combigrid_convenience.patch
./bin/spack compiler find

# if there is an intel compiler on the system, use it
if ./bin/spack spec sgpp@adaptive_combigrid_convenience+combigrid+python%intel^swig@3.0.12; then
    # in case of version conflicts, may need to restrict dependency versions further
    # e.g., on my laptop:
    #./bin/spack spec sgpp@adaptive_combigrid_convenience+combigrid+python^py-numpy@1.17:1.18^py-scipy@1.3:1.5^py-setuptools@:40

    # at this point, you can tell spack which components are already installed by
    # editing ./etc/spack/packages.yaml -- this avoids re-installing everything
    # cf. https://spack.readthedocs.io/en/latest/build_settings.html
    # this here may be useful too:
    # ./bin/spack external find
    # if you are on cobra, run this to find python packages:
    # git apply ../cobra_spack_python.patch

    # this here can take a long time:
    ./bin/spack install sgpp@adaptive_combigrid_convenience+combigrid+python%intel^swig@3.0.12
else
    ./bin/spack spec sgpp@adaptive_combigrid_convenience+combigrid+python^swig@3.0.12 &&
        ./bin/spack install sgpp@adaptive_combigrid_convenience+combigrid+python^swig@3.0.12
fi

# load sgpp paths -- do this every time you need pysgpp on your PATH
. ./share/spack/setup-env.sh
spack load sgpp@adaptive_combigrid_convenience

cd ..
git clone git@gitlab.mpcdf.mpg.de:tpoll/gene-diag.git
