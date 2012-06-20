MultiDimEvol
============

*Simulating the multi-dimensional adaptive evolution of a feed-forward network*

This repository contains the code and analyses relating to:

Bullaughey, K. (2012) Multidimensional adaptive evolution of a feed-forward network and the illusion of compensation. Evolution (in press)

Dependencies & Runtime Environment
----------------------------------

I have only tested this code on OS X 10.6 and Redhat Enterprise Linux 5.4. 

The code base is a combination of R, C++, and Bash scripts. The C++ portion is compiled into a library that can be dynamically loaded into R. 

The R code depends on the following libraries: `grid`, `xtable`

Installation
------------

To compile the R library:

    cd src/c
    make

This will produce a file `reg_xyz.so`, which can be loaded into an R session with the command:

    dyn.load("reg_xyz.so")

Compiling the manuscript
------------------------

I have also included the source for the manuscript, as submitted to Evolution. Both the main text and supplement are written in LaTeX. Compiling it will require `pdflatex`, `bibtex`, `dvips` `gs` (ghostscript) and `make`.  

Simply go into the manuscript directory and run make:

    cd manuscript
    make

This will produce two PDFs:

    build/bullaughey-2012_evolution.pdf
    build/bullaughey-2012_evolution-supporting_information.pdf

Here are the submitted versions of the [main text][mt] and [supplement][sup].

[mt]: https://github.com/kbullaughey/multidimevol/blob/master/manuscript/submitted/bullaughey-2012_evolution.pdf?raw=true
[sup]: https://github.com/kbullaughey/multidimevol/blob/master/manuscript/submitted/bullaughey-2012_evolution-supporting_information.pdf?raw=true

Data behind the figures and tables
----------------------------------

**Table S1**

