# Part 2: Introduction
For our first practical session, we will be using `BEAST` to reconstruct the evolutionary 
dynamics of influenza. First, we will need to convert the input data, which is NEXUS format, 
into BEAST XML format. For that purpose, we will use `BEAUti`. Both software are available for 
Mac OS X, Windows, and UNIX/Linux OSs when you download the latest stable version 
of `BEAST V1.10.4` (see [here](https://github.com/beast-dev/beast-mcmc/releases/tag/v1.10.4)). 

To visualise the output file of a Bayesian MCMC software (for this second part we will be 
using `BEAST` and `MCMCtree`), we can use [`Tracer`](https://github.com/beast-dev/tracer/). 
This program is very useful to summarise, both visually and quantitatively, the samples 
that have been collected during the MCMC. It can also be used as a diagnostic tool as 
you can plot the traces of the values sampled for the parameters of interest, it computes 
the ESS for each parameter, etc. 

Last, we wil use [`FigTree`](https://github.com/rambaut/figtree/releases)
to display the phylogenies that we will be working with during this tutorial.
