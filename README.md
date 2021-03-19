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

## Running the GUI

1. Choose a directory you want to put the code `/path/to/root/`

2. Choose a directory where you want to keep data `/path/to/data/`

3. Go to the code directory using  cd `/path/to/root/`

4. Clone the git repo
```
git clone git@github.com:BenSmithers/Cascade-Smearing.git
```
5. Rename the downloaded folder, `mv Cascade-Smearing cascade`. Without doing this, the code will be looking for a module of the wrong name. 

6. Go into the newly cloned repository, Cascade-Smearing, and edit the datapath line in `config.json` to read
  ```
  "datapath": "/path/to/data",
  ```
replacing `/path/to/data` with your actual path. This tells my code where to look for the data.

7. Edit your .bashrc, adding in the line
```
export PYTHONPATH=$PYTHONPATH:/path/to/root/
```
This tells python how to find the code.

8. Now you need to download the data. If you have access to the cobalts, go to /path/to/data and run
```
scp username@cobalt.icecube.wisc.edu:/data/user/bsmithers/data/cascades/hg_sib/flux_data_* .
scp username@cobalt.icecube.wisc.edu:/data/user/bsmithers/data/cascades/hg_sib/downsized_* .
```
Note: this might take a while! I have a few other folders for fluxes from different MCEq fluxes. 

9. You should now be able to run the gui, after restarting whatever terminal you're using, by going to `/path/to/root/cascade/plotting` and running 
```
python3 explo.py
```

Note: I haven't tested this in Python2 and have no intentions of supporting Python2. 
I'm pretty sure PyQt5 and some class decorators work differently. 
If you *are*  still using Python2... *why?*

### Updates

 - There is new a drop-down menu to turn different neutrino fluxes on and off in the displayed ratios. 
 
 - There is a new slider to adjust the colorbar scale. It will always be centered on 1.0, but can have a width anywhere from (0.99-1.01) up to (0.00-2.00)
 
 - You can switch between reconstructed and raw fluxes (with cross-sections) using the checkbox at the bottom of the screen. Switching to raw fluxes takes a long time right now, and so the sliders are disabled when in raw flux mode. 

## Links

These aren't realated, I'm just trying to keep track of these links for my own sake. 

https://github.com/IceCubeOpenSource/HESE_data_release

https://github.com/fatlamb/verosimilitud


