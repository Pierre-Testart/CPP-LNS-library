## Logarithmic Number System library

This project was carried out by Pierre Testart, student at Politecnico di Milano.  
Its aim was to create a header-only library to support the Logarithmic Number System in C++, and then to integrate benchmarks (from [Axbench](http://axbench.org/)) to measure accuracy and performance.  
This is the header-only library, the benchmarks will be uploaded in a different Github repository.  
C++14 is the required standard to use the library.

### How to use the LNS types

This library provides a class `lns_t<I, F, A>`, inside the `lns` namespace.  
The 3 template parameters are integers.
* I is the number of integer bits used for the exponent of the number (including sign bit). It must be at least 1.
* F is the number of fractional bits used for the exponent of the number. It must be positive.
* A is the level of approximation of the lookup tables used for faster LNS addition. With a greater value of A, lookup tables will be lighter in memory, but less accurate. The precision of tables is such that the maximum difference between 2 consecutive values is no more than 2^A. If A is set to -1, faster addition is disabled, and no lookup tables are used. There is a default value for A based on the value of F, making this parameter optional.

Since the exponent is an integer, there is a constraint that I+F must be either 16, 32 or 64.  
There are predefined types for each exponent size: `lns16_t`, `lns32_t` and `lns64_t`, defined in the `lns` namespace.  
Note that the size of a LNS variable in memory is not the size of the exponent, there is an additional byte used for storing sign and flags for special values.
