# Cascade Smearing

## History

This little project was initially started in my `Analysis` repository, and so for a full history of the commit history you can check out [this link](https://github.com/BenSmithers/Analysis/tree/master/nusquids_stuff/convolve). 
The project eventually evolved beyond the point where it was still reasonable to keep it as a sub-folder in a bigger repo. 

## Introduction

What this does is it uses MCEq to get neutrino fluxes at IceCube. 
Then, we use differential fluxes from NUSQuIDS to get differential fluxes for energy deposited vs true energy vs neutrino zenith angle. 
These then are combined with reconstruction data from the 7yrs Cascades paper to get a 4D array of fluxes:
    - True Event Energy
    - Reconstructed Energy
    - True cos(zenith) angle of event
    - Reconstructed cos(Zenith) angle of event

Then, we can take a tuple of reconstructed event params (energy, cos(Zenith)) and do a binlinear interpolation of the 2D slice (true E, true cz) around solved points in that 4D lattice. That yields a 2D array of likely truth values: their relative values indicating relative likelihoods of getting some truth value given reconstructed ones. 

## The files:

What do each of the files do, at a glance? 
    - convolve.cpp : this calculates the fluxes from mceq, and propagates them through the Earth using nusquids
    - cross_section_test : this just holds some macros for accessing nusquids cross section functions
    - deporeco : handles the probabilities for reconstructing an energy given a deposited energy 
    - deposit : does the work for taking the neutrino fluxes and building that true E vs. depo E array 
    - imagescan : not really used anymore. Scans a picture from the cascades paper to get probabilities for reconstructing some energy given a deposited energy 
    - mceq_flux : the mceq fluxes
    - nus_compile.sh : this is a bash script for compiling convolve.cpp
    - nus_utils : more nusquids utilities 
    - plot_contour: plots the 2D likelihood plots discussed in the Introduction
    - plot_fluxes: plots some of the outputs from `deposit.py`, used now mainly as tests 
    - reconstruct: makes the 4D lattice discussed in the Inroduction 
    - tau_funcs: focuses on the weirder tau-energy-deposition probabilities 
    - utils: a bunch of cool utility functions! 
