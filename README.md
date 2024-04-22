# MSci Project (in Python)


A group of Python scripts, which encapsulate different approximations made for modelling the formation of giant planets in protoplanetary discs.


## Data

The data the scripts use.


## Plots

The outputted plots from the scripts of this project.


## Scripts

### IsoT_Model.py
A class which allows a user to input initial conditions for the instantiation of a giant planet with an isothermal envelope. This simple model allowed us to test whether our code was working as expected as we can compare numerical outputs to their analytic equivalents.
### IsoD_Model.py
The same as IsoT_Model.py, but for an envelope with constant density.
### IsoO_Model.py
A class similar to above, but only uses a constant opacity assumption for the envelope. There are methods allowing users to run the equations governing the envelopes structure (Planet_Formation()), calculate the convergence stability of the code (L2norm()), and calculate the parameters which best describe the density profile (DensityFit()) - this would be used when determining the ideal initial conditions for the grid of mass values such that code stability is ensured.
### L_Finder_IsoO.py
A function which allows a user to find the luminosity value which causes the inner envelope radius to fall on a given giant planetary core. This luminosity can be used to determine how much mass the planet can accrete per unit time (thus evolve it).
### Evolver_IsoO.py
A function which allows a user to evolve the envelope of a giant planet by increasing its mass according to its luminosity value and rerunning the IsoO_Model method “Planet_Formation()”, and so on. This function can thus provide us with an mass accretion rate.
### GPF_Movshovitz_Opacity_Forced.py
This is the same as “IsoO_Model.py”, but uses an opacity profile taken from the Movshovitz et al. 2010 paper. The “_Forced” refers to the fact that the convective/radiative boundaries are hard-coded in.
### GPF_Movshovitz_Opacity.py
The same as “GPF_Movshovitz_Opacity_Forced.py” but where the convective/radiative boundaries arise naturally from the Schwarzschild criterion. This tested our implementation of the Schwarzschild criterion was accurate is our plots could be compared to the Movshovitz 2010 paper plots.
### Freedman_Opacities_Function.py
A function to return the opacities calculated from the analytic fit from the Freedman et al 2014 paper.
### GPF_Freedman_Opacity.py
A class much like “GPF_Movshovitz_Opacity.py”, but using the Freedman opacity function to describe the opacity. There now is included a Brent’s method to find the temperature and opacity value for a given grid point simulataneously (as the opacity and temperature are now dependent upon one another).
### L_Finder_FreO.py
The luminosity finder for the Freedman opacity approximation.
### Evolver_FreO.py
The evolver for the Freedman opacity approximation.
### Combined_Opacities_Function.py
A script of functions which accumulates into the last function (Combined_Opacity()) which allows for a user to calculate an opacity from either the dust or gas depending on what temperature region the opacity is being calculated for (i.e. whether the dust has sublimed). The dust opacities are taken from Mordisini et al. 2012 paper b.
### GPF_Combined_Opacity.py
A class much like “GPF_Freedman_Opacity.py”, but using the “Combined” opacity function to describe the opacity.
### L_Finder_ComO.py
The luminosity finder for the “Combined” opacity approximation.
### Evolver_ComO.py
The evolver for the “Combined” opacity approximation.
### Maximum_Accretion_Rate_Calculator.py
A function which calculates the maximum rate at which matter can be accreted by a forming giant planet, dependent upon how fast the surrounding disc can provide material. Different assumptions were made, informed from Tanigawa et al. 2012 and Lissauer et al. 2009, with the surface density profile of the disc informed by Crida et al. 2006.
### Data_Loader.py
Functions that extract and format data from external files into a more easily readable format for other scripts.
### Mass_Grid_Solver.py
Functions that take the parameters output from DensityFit() methods in aforementioned classes and outputs a grid of mass values which would provide a logarithmically spaced grid of radius values for each mass value. (This would be repeatedly found until the mass grid doesn’t change - most stable grid using this method).
### Convergence_Tester.py
A script with many functions to test the convergent properties of the different numerical integrations throughout this project.
### Tau_Line_Finder.py
A function to find where the optical depth is equal to 1 for a given envelope (used for plots).
### Plotter.py
A list of functions allowing a user to plot data outputted from the previous scripts.
### Test.py
The script used to run all the previous classes and functions to output the plots in the “Plots” folder.
