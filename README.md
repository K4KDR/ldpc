# LDPC
*This project is not yet finished. Especially the decoder is known to be buggy, 
only use this for testing*

C++ project to handle LDPC messages in the MOVE-II project.

This repository contains:
1. The core library
2. Unittests/ example applications
3. (Optional) Application to compute systematic generator matrices

Each part is installed by cmake. The entire project can be build at once with
````
cd ldpc
mkdir build && cd build
cmake ../
make make test
sudo make install
````

In order to include the generator matrix computation programm (which requires
NTL), change the cmake command in the above procedure to
`cmake -DBUILD_GENERATOR=On ../`.

## Core library (libldpc) ##
The core library is a shared library that exports the two classes encoder and
decoder. They can be included with `#include <ldpc/encoder.h>` and
`#include <ldpc/decoder.h>`.

## Example applications
The unittests in the tests folder serve as demonstrations how to use the
library. They can be run from the build folder with `make test`.

## Application to compute systematic generator matrix
The application `ldpc_compute_generator` computes a generator matrix from a
given parity check matrix (in alist format). The application assumes that the
code is a systematic code, where the input bits are send before the computed
check bits.

This application depends on the NTL library, that has to be installed.

In order to generate the binary output file (*.gen) that is required by the ldpc
library, the binary needs to be called with 3 arguments. The first argument is
the path to the alist file defining the parity check matrix. The second argument
is a path where an ASCII version of the computed generator matrix is written to.
The third argument is the path where the *.gen file will be written to.

````
ldpc_compute_generator paritycheck_matrix.a generator.txt generator.gen
````
