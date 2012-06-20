MultiDimEvol
============

*Simulating the multi-dimensional adaptive evolution of a feed-forward network*

This repository contains the code and analyses relating to:

Bullaughey, K. (2012) Multidimensional adaptive evolution of a feed-forward network and the illusion of compensation. Evolution (in press)

Dependencies & Runtime Environment
==================================

I have only tested this code on OS X 10.6 and Redhat Enterprise Linux 5.4. 

The code base is a combination of R, C++, and Bash scripts. The C++ portion is compiled into a library that can be dynamically loaded into R. 

Installation
============

To compile the R library:

    cd src/c
    make

This will produce a file `reg_xyz.so`
