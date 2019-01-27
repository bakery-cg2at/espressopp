This is a modified version of ESPResSo++ to work with reverse mapping method that can be found here:
https://github.com/bakery-cg2at/bakery


We closely based on the upstream version of ESPResSo++, at the level of API and the implementation of the
potential.
The additional part lays in the AdResS section, which was reimplemented from scratch. If you use in your research any part of the new code (namely ``FixedVSList``, ``VelocityVerletHybrid``, ``VerlelListHybrid``), then we kindly ask to cite the following paper:

```
@article{doi:10.1021/acs.jctc.6b00595,
  author = {Krajniak, Jakub and Pandiyan, Sudharsan and Nies, Eric and Samaey, Giovanni}
  title = {Generic Adaptive Resolution Method for Reverse Mapping of Polymers from Coarse-Grained to Atomistic Descriptions},
  journal = {Journal of Chemical Theory and Computation},
  doi = {10.1021/acs.jctc.6b00595},
  note ={PMID: 27685340},
  URL = { http://dx.doi.org/10.1021/acs.jctc.6b00595},
}
```

Please note, that according to the software licence, if your published work uses any modification of this package, you are obliged to make the modification publicly available. We are happy to accept any of the improvements to the code base.


# Install

In order to install the package, please take the following steps

```
# cd espressopp
# mkdir build
# cd build
# cmake ..
# make
# source ESPRC
```

this will compile the package and sourcing `ESPRC` will setup all nessecary envrionmental variables.
You can check the correctness of the installed package by importing `esressopp` module:
```bash
python -c "import espressopp"
```
if the import does not return any error, then the Python package is set up correctly.
For further checkup, you can run unit-tests: `make test`.

## Install on HPC with Intel compiler

In order to use Intel compiler toolchain you have to explicitly put the path to the `icc` and `icpc` binaries in the `cmake` command:

`cmake -DCMAKE_CXX_COMPILER=<path to icpc> -DCMAKE_C_COMPILER=<path to icc> ..`

**Important note**

Usually, HPC provides multiple toolchains, e.g. for gcc or intel compilers. This is mostly managed by some sort of [software environment modules](http://modules.sourceforge.net/).
Keep in mind, that if you use intel toolchain then all of the libraries have to be compiled with the same compiler, e.g. boost, mpi or fftw3 libraries have to 
compiled with intel compiler.



## Troubleshooting

### Could NOT find FFTW3 (missing: FFTW3_LIBRARIES FFTW3_INCLUDES)

CMake could not find FFTW3 library in any of the standard system location. You have to provide explicitly path to the FFTW3 library:

`cmake -DFFTW3_LIBRARIES=<path to libfftw3.so> -DFFTW3_INCLUDES=<path to include>`
e.g.
`cmake -DFFTW3_LIBRARIES=/usr/lib/libfftw3.so -DFFTW3_INCLUDS=/usr/include`


# Issues

Report bugs on the [github issues site](https://github.com/bakery-cg2at/espressopp/issues)
