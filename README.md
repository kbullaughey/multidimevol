MultiDimEvol
============

*Simulating the multi-dimensional adaptive evolution of a feed-forward network*

This repository contains the code and analyses relating to:

Bullaughey, K. (2012) Multidimensional adaptive evolution of a feed-forward network and the illusion of compensation. Evolution (in press)

The authorative copy of this source can be found at: https://github.com/kbullaughey/multidimevol

Getting the code
----------------

The preferred way to download the source, data, and plots, is to clone the git repository like so:

    git clone git://github.com/kbullaughey/multidimevol.git

This will create a directory `multidimevol` containing everything. 

**Please note that this repository is rather large (hundreds of MB), because it includes all the simulated data.** 

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

Running the analyses
--------------------

**Note:** Unfortunately, I did not always set a seed for the random number generators, and so the results in the paper are not exactly reproducable from the code contained in this repository. Nonetheless, I do not expect that any of the results are outliers, and so the running the analyses contained here should produce essentially the same results.

All of the analyses scripts are contained in the directory `analyses`. It's easiest if these scripts are run from within this directory because they use a relative path to locate the `src` directory.

In general, I pass a parameter `run` to each script. This can be any string that can appear in a file name, and I use it to version my analyses. Throughout the analyses, the run corresponding to the submitted version is labeled `2011_10_20`. All the data files and plots from this run are in the git repository.

**Producing Figure 4 and Table S1**

To make a new version of Figure 4 (with a new run name), there are two steps. The first performs the simulations:

    R --vanilla --args --run=2012_06_20 < random_starts-run.r 

And the second step, produces Figure 4 as well as the data for Table S1:

    R --vanilla --args --run=2012_06_20 < random_starts-curves-grid.r

**Producing Figure 5**

To make Figure 5, run the following:

    R --vanilla --args --run=2012_06_20 < two_environments.r

This will produce a PDF:

    analyses/plots/two_environments-2012_06_20.pdf

And several configs:

    analyses/configs/two_environments-start-2012_06_20.rconf
    analyses/configs/two_environments-low_noise_optimum-2012_06_20.rconf
    analyses/configs/two_environments-high_noise_optimum-2012_06_20.rconf

Corresponding to a random starting point in the low-noise environment (start), the optimum when the network is allowed to evolve adaptive in the low-noise environment (low_noise_optimum) and the optimum after the low-noise optimal network is switched into the high-noise environment and allowed to evolve adaptively until it reaches a maximum.

**Producing Figure 6**

Figure 6 involved running 100 replicate simulations, staring from a particular feed-forward network that had reached a maximum by evolving adaptively to the low-noise environment, and evolving it in the high-noise environment under the assumptions of discrete mutaitonal effects, each of which only affects a single trait.

The initial feed-forward network that is adapted to the low-noise environment is given in:

    configs/discrete_effect_start.rconf

And the simulation paremters I used are given in:

    configs/discrete_effects.sconf

Given I had a cluster of compute nodes available running Sun Grid Engine, I submitted the jobs to the cluster:

    for i in `seq 1 100`; do 
      qsub discrete_effects-submit.sh 2012_06_20 configs/discrete_effect_start.rconf configs/discrete_effects.sconf
    done

These jobs are not particularly computationally intensive (taking less than 1 minute each on contemporary hardware), and so they could be run sequentially. However, the shell script assumes the existence of an environment variable, `JOB_ID`, which one would need to set uniquely for each repliate, so the files do not overwrite each other.

The above jobs create a series of directories containing the output of the following form:

    analyses/data/discrete_effects-simulations/2012_06_20/$JOB_ID

Each run is plotted in a PDF named `2012_06_20-$JOB_ID-selection.pdf` in the corresponding directory. Here is [one example for job 6798066][selplot].

[selplot]: https://github.com/kbullaughey/multidimevol/blob/master/analyses/data/discrete_effects-simulations/2012_06_20/6798066/2012_06_20-6798066-selection.pdf?raw=true

