# EAGLES
## Estimating AGes from Lithium Equivalent widthS

This is the software described in the paper "Estimating Ages from Lithium Equivalent Widths (EAGLES)" by Jeffries et al. (2023), MNRAS submitted. 
The code implements an empirical model that predicts the lithium equivalent width (EW) of a star (from the Li I 6708A line), as a function of its age and effective temperature ($T_{eff}$). This was developed by fitting to a training dataset consisting of around 6000 stars in 53 open clusters that were observed as part of the Gaia-ESO survey ([Gilmore et al. 2022](https://ui.adsabs.harvard.edu/abs/2022A%26A...666A.120G/abstract); [Randich et al. 2022](https://ui.adsabs.harvard.edu/abs/2022A%26A...666A.121R/abstract); [Jackson et al. 2022](https://ui.adsabs.harvard.edu/abs/2022MNRAS.509.1664J/abstract)) and is applicable to stars with $3000 < T_{eff} < 6500$ K and $-0.3 < [Fe/H]  < 0.2$.

This model is used to compute the age probability distribution for a star with a given EW and Teff, subject to an age probability prior that may be flat in age or flat in log age. If the data for more than one star is entered, then the option exists to treat this as a cluster and determine the age probability distribution for the ensemble. The code will produce estimates of the most probable age, uncertainties and the median age, output files consisting of probability plots, best-fit isochrone plots and tables of the posterior age probability distribution(s).


## Running the code

The code consists of a simple Python 3 script (eagles.py), which contains an extensive help file that details the data input format and the command line options. ThIs can be viewed at the head of the code or by executing it as "**<script> -h**"
  
**input.dat** is an example input data file containing the data for three stars which might be part of a cluster
  
To test whether the code is working use  
  **"<script> input.dat output -c -s"**,  
  which should report the following  
  
 
![image](https://user-images.githubusercontent.com/104770145/218447810-62c683af-cb86-4c0a-82e6-b1f58c523ebb.png)


 Cluster of 3 stars  
 chi-squared of fit =   0.46  
 most probable log (Age/yr) = 7.277 +0.135/-0.630  
 most probable Age (Myr) =   18.9 +   6.9/-  14.5  
 median log (Age/yr) = 7.112  
 median Age (Myr) =   13.0  
  
and produce the output files
  
*  output_prob.csv       - plot of the combined age probability distribution
*  output_iso.pdf        - plot of the data and best-fit isochrone in the EW vs Teff plane
*  output_pos.csv        - combined posterior probability distribution for the dataset
*  output.csv            - summary results for the three stars and result for the cluster
  
  
 ## Additional Scripts/Files
 
 * eagles_iso.py   - script to produce model isochrones of EWLi vs Teff, plot and save them as ascii files
 * eagles_iso.zip  - a zip file containing isochrones at 5, 10, 15, 20, 30, 50, 100, 200, 300, 500, 1000, 2000, 3000, 5000 Myr. Any other ages can be produced using eagles_iso.py
 
 ![image](https://user-images.githubusercontent.com/104770145/218448006-21132158-63b4-49a3-a6f1-e8cfcc28de37.png)

 
  # Contact
  To register an interest, request clarifications or report bugs - email r.d.jeffries@keele.ac.uk
