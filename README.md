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

The R code depends on the following libraries: `grid`, `xtable`.

Other dependencies include `perl` and basic unix commands like `ls`.

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

A note about compute clusters
-----------------------------

When running these analyses, I had access to a large compute cluster running Sun Grid Engine. Therefore the workflow I describe involves launching jobs using the `qsub` command. Each job has a submission script written in `bash` that sets up a memory allocation and output files for `stdin` and `stderr`. These submission scripts assume the existence of an environment variable, `JOB_ID`, which one would need to set uniquely for each repliate, so the files do not overwrite each other. 

For folks not running on a cluster, I have provided a mock `qsub` command, `src/pl/qsub`, which simply picks a random number for the `JOB_ID` and then runs the submission script directly in a shell environment with this variable defined. 

Simply copy `src/pl/qsub` into a directory that is in your path, make it executable (if it's not already), and you should be able to run the workflows below, albeit slower than when executed in parallel.

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

These jobs are not particularly computationally intensive (taking less than 1 minute each on contemporary hardware), and so they could be run sequentially. See the above section on compute clusters. 

The above jobs create a series of directories containing the output of the following form:

    analyses/data/discrete_effects-simulations/2012_06_20/$JOB_ID

Each run is plotted in a PDF named `2012_06_20-$JOB_ID-selection.pdf` in the corresponding directory. Here is [one example for job 6798066][selplot].

[selplot]: https://github.com/kbullaughey/multidimevol/blob/master/analyses/data/discrete_effects-simulations/2012_06_20/6798066/2012_06_20-6798066-selection.pdf?raw=true

**Producing Figure 7**

The analysis presented in Figure 7 involves all non-empty subsets of the five evolvable traits. Here I perform some preliminary setup, determining all non-empty subsets, and producing a shell script to launch the jobs:

    R --vanilla --args --run=2012_06_20 < dimension_analysis-setup.r

The shell script produced can be run as follows:

    bash dimension_analysis-launch-2012_06_20.sh

This will create `.rimage` files with the naming convention:

    analyses/data/dimension_analysis/2012_06_20/dimsim-<params>.rimage

The output requires some post-processing:

    ./dimension_analysis-post_process.sh 2012_06_20

This will create `.rimage` files with the naming convention:

    analyses/data/dimension_analysis/2012_06_20/dimsim-<params>-biggest_reversals.rimage

At this point, Figure 7 can be generated as follows:

    R --vanilla --args --run=2012_06_20 < dimension_analysis-plot.r

Which will produce a PDF at this location:

    analyses/plots/dimension_analysis-2012_06_20.pdf

Which looks like [this][da].

[da]: https://github.com/kbullaughey/multidimevol/blob/master/analyses/plots/dimension_analysis-2012_06_20.pdf

** Producing Figure S1**

After the dimension analysis (above) is run, Figure S1 can be produced as follows:

    R --vanilla --args --run=2012_06_20 < dimension_analysis-sizes_by_dimension.r

This results in a pdf at this location:

    analysis/plots/reversal_sizes_by_dimension-2012_06_20.pdf

