# LDPC library
*This library is not yet finished. Especially the decoder is known to be buggy, 
only use this for testing*

C++ library to handle ldpc messages.

This repository contains:
1. The core library
2. Example applications
3. Sample application to compute systematic generator matrices

## Core library
The library build process is controlled by make. A typical setup can be performed with
````
cd libldpc
make
sudo make install
````

## Example applications
Example applications using the library are available in the tests directory. They 
can be build with the provided Makefile, placing the resulting binaries in the tests folder as well.
````
cd libldpc/tests
make
./test_encoder
./test_decoder
./test_chain
````

## Application to compute systematic generator matrix
Based on a parity check matrix provided in the alist format the programm computes 
a generator matrix, assuming that the code is a systematic code, where the input 
bits are send before the computed check bits.

This application depends on the NTL library, that has to be installed.

The build process is controlled by Make.
````
cd libldpc/compute_gen
make all
sudo make install
````

In order to generate the binary output file (*.gen) that is required by the ldpc 
library, the ldpc_compute_gen binary needs to be called with 3 arguments. The 
first argument is the path to the alist file defining the parity check matrix. 
The second argument is a path where an ASCII version of the computed generator 
matrix is written to. The third argument is the path where the *.gen file will be 
written to.
````
ldpc_compute_gen paritycheck_matrix.a generator.txt generator.gen
````
