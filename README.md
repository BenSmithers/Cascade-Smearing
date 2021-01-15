# Cascade Smearing

## History

This little project was initially started in my `Analysis` repository, and so for a full history of the commit history you can check out [this link](https://github.com/BenSmithers/Analysis/tree/master/nusquids_stuff/convolve). 
The project eventually evolved beyond the point where it was still reasonable to keep it as a sub-folder in a bigger repo. 

## Introduction

What this does is it uses MCEq to get neutrino fluxes at IceCube given a set of sterile oscillation parameters. 
Then, we use differential fluxes from NUSQuIDS to get differential fluxes for energy deposited vs true energy vs neutrino zenith angle. 
These then are combined with reconstruction data from the 7yrs Cascades paper to get a 4D array of fluxes:
    - True Event Energy
    - Reconstructed Energy
    - True cos(zenith) angle of event
    - Reconstructed cos(Zenith) angle of event

Then, we can take a tuple of reconstructed event params (energy, cos(Zenith)) and do a binlinear interpolation of the 2D slice (true E, true cz) around solved points in that 4D lattice. That yields a 2D array of likely truth values: their relative values indicating relative likelihoods of getting some truth value given reconstructed ones.

We also make just do a straight up convolution to get fluxes in reconstruction space. This is what we really like.  

## The files:

What do each of the files do, at a glance?

    - config.json : a JSON-formatted text file providing configuration for the run  


    - cross_section_test : this just holds some macros for accessing nusquids cross section functions.
It also keeps the cross sections in an object

    - deporeco : handles the probabilities for reconstructing an energy given a deposited energy 

    - deposit : does the work for taking the neutrino fluxes and building that true E vs. depo E array 

    - nus_utils : more nusquids utilities

    - raw_fluxes : generates the true neutrino fluxes at the atmosphere and at the detector using MCEq and nuSQuIDS  

    - reconstruct : makes the 4D lattice discussed in the Inroduction

    - sensitivity : right now this just makes some plots, more to come soon 

    - tau_funcs: focuses on the weirder tau-energy-deposition probabilities 

    - utils: a bunch of cool utility functions! 



## Links

https://github.com/IceCubeOpenSource/HESE_data_release

https://github.com/fatlamb/verosimilitud


