# SRIM_Width_finder
Collection of python codes used to analyse multiple SRIM files to find the correct thickness

By iteration across multiple distances, a collection of particles transmitted through a foil with an energy less than 5 keV can be made.
A measurement of these across multiple thickness returns a gaussian distribtion, which we use to find the "Optimal" thickness. 
Bootstrap resampling is used to calculate percentage errors at given thickness'.

Srim_Distance calculates the gaussian curve fit across multiple data sets for sub 5 keV transmission, and produces plots that visually show this.
Srim_Plotter plots the optimised foil distributions, along with the errors resulting from bootstrap resampling.

gaussians contains a collection of useful functions which are used within these codes.

Use Srim_Distance first to find the optimised thickness, then run a simulation at that thickness and adjust the code accordingly to reflect the optimised thickness.

Afterwards use Srim_Plotter to recieve useful plots that demonstrate the efficiency of each foil

