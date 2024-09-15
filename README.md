# EAGLES
## Estimating AGes from Lithium Equivalent widthS

This is the software described in the paper "Estimating Ages from Lithium Equivalent Widths (EAGLES)" by Jeffries et al. (2023), MNRAS, 523, 802 ([https://arxiv.org/abs/2304.12197](https://arxiv.org/abs/2304.12197)). 
The code implements an empirical model that predicts the lithium equivalent width (EW) of a star (from the Li I 6708A line), as a function of its age and effective temperature ($T_{eff}$). This was developed by fitting to a training dataset consisting of around 6000 stars in 52 open clusters that were observed as part of the Gaia-ESO survey ([Gilmore et al. 2022](https://ui.adsabs.harvard.edu/abs/2022A%26A...666A.120G/abstract); [Randich et al. 2022](https://ui.adsabs.harvard.edu/abs/2022A%26A...666A.121R/abstract); [Jackson et al. 2022](https://ui.adsabs.harvard.edu/abs/2022MNRAS.509.1664J/abstract)) and is applicable to stars with $3000 < T_{eff} < 6500$ K and $-0.3 <$ [Fe/H]$  < 0.2$.

This model is used to compute the age probability distribution for a star with a given EW and Teff, subject to an age probability prior that may be flat in age or flat in log age. If the data for more than one star are entered, then the option exists to treat this as a cluster and determine the age probability distribution for the ensemble. The code will produce estimates of the most probable age, uncertainties and the median age, output files consisting of probability plots, best-fit isochrone plots and tables of the posterior age probability distribution(s).

A new version 2.0 of the code implements a machine learning model for the relationship between EWLi, $T_{eff}$ and age (see  [Weaver, Jeffries and Jackson 2024, arXiv 2409.07523](https://arxiv.org/abs/2409.07523)
and below).

## Running the code

The code consists of a simple Python 3 script (eagles.py), which contains an extensive help file that details the data input format and the command line options. ThIs can be viewed at the head of the code or by executing it as "**<script> -h**"
  
**input.dat** is an example input data file containing the data for three stars which might be part of a cluster
  
To test whether the code is working use  
  **"<script> input.dat output -c -s"**,  
  which should report the following  
  
![github1](https://user-images.githubusercontent.com/104770145/234336618-d72f732a-8d97-4878-b6b8-87b2cb07b69a.png)



 Cluster of 3 stars  
 chi-squared of fit =   0.46  
 most probable log (Age/yr) = 7.272 +0.135/-0.630  
 most probable Age (Myr) =   18.7 +   6.8/-  14.3  
 median log (Age/yr) = 7.112  
 median Age (Myr) =   13.0  
  
and produce the output files
  
*  output_prob.pdf       - plot of the combined age probability distribution
*  output_iso.pdf        - plot of the data and best-fit isochrone in the EW vs Teff plane
*  output_pos.csv        - combined posterior probability distribution for the dataset
*  output.csv            - summary results for the three stars and result for the cluster
  
  
 ## Additional Scripts/Files
 
 * eagles_iso.py   - script to produce model isochrones of EWLi vs Teff, plot and save them as ascii files
 * eagles_iso.zip  - a zip file containing isochrones at 5, 10, 15, 20, 30, 50, 100, 200, 300, 500, 1000, 2000, 3000, 5000 Myr. Any other ages can be produced using eagles_iso.py
 
 ![image](https://user-images.githubusercontent.com/104770145/218448006-21132158-63b4-49a3-a6f1-e8cfcc28de37.png)


## EAGLES V2

A newer version of EAGLES is included in the zip file 'eaglesv2_0.zip'. This new version inclues an option for an artificial neural network (ANN) model of the relationship between EWLi (and its intrinsic dispersion), $T_{eff}$ and (log) age ([Weaver, Jeffries and Jackson 2024, arXiv 2409.07523](https://arxiv.org/abs/2409.07523)). The zip file contains a new script 'eaglesv2_0.py' along with a pre-computed grid of EWLi (and its intrinsic dispersion) as a function of $T_{eff}$ and age, a readme file containing a description of the changes in version 2.0 and the same 'input.dat' test file.

The machine learning model provides better accuracy in reproducing the relationship between EWLi and its dispersion with age, free from the constraints of an arbitrary analytical model (see Weaver et al. 2024 for a full discussion).

To run the script using the ANN model, simply use the -m option on the command line, making sure that the pre-computed grid is in the same folder as the script. If the -m flag is omitted then the model used is the same analytic model described by Jeffries et al. (2023) and implemented in version 1.0.

Running on the example 'input.dat' file

 **"<script> input.dat output -m -c -s"**,  
  which should output the following  

![eaglesv2](https://github.com/robdjeff/eagles/assets/104770145/31fe7a94-fe59-4b09-8961-aff1dbf39a4f)


 **EAGLES V2.0**  

 No additional Teff errors in input file  

 Using ANN model:  

 Setting log(age) limits to 6.0 - 10.1

 Cluster of 3 stars  
 chi-squared of fit =   0.20  
 most probable log (Age/yr) = 7.315 +0.160/-0.640  
 most probable Age (Myr) =   20.7 +   9.2/-  15.9  
 median log (Age/yr) = 7.095  
 median Age (Myr) =   12.4  

and produce the output files  
  
*  output_prob.pdf       - plot of the combined age probability distribution
*  output_iso.pdf        - plot of the data and best-fit isochrone in the EW vs Teff plane
*  output_pos.csv        - combined posterior probability distribution for the dataset
*  output.csv            - summary results for the three stars and result for the cluster

  # Citation
  If you make use of the EAGLES code(s) we would appreciate that you cite [Jeffries et al. 2023, MNRAS, 523, 802](https://academic.oup.com/mnras/article/523/1/802/7147327).

  If you use Version 2 od the code (the neural network model), you should also cite [Weaver, Jeffries and Jackson 2024, arXiv 2409.07523](https://arxiv.org/abs/2409.07523).
  # Contact
  To register an interest, request clarifications or report bugs - email r.d.jeffries@keele.ac.uk
