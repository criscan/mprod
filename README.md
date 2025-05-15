# mprod (GSPM)
Generalized Surplus Production Model (GSPM)
(by Cristian M. Canales)

Out=mprod_fits(name, graf_opt, tab_opt)

name     = "name_file.xlsx"
graf_opt = T or F (to produce model fit and population variables graphics)
tab_opt  = T or F (to produce output tables)


An Excel file is read with data separated into two spreadsheets: one with catches and CPUE by year, and the other with parameters priors. The model estimates provide different interest variables and parameters. Both files must be at same location folder

The second spreadsheet considers the start parameters. The initial biomass as carrying capacity (K) could be established as K=0.8 * Sum of historical landings, r is the growth population rate (depending on the species biology), sigma is the standard deviation of data vs model, rho is the asymmetry rate of growth function (1=symmetric Shaefer model), and p0 a guess on population depletion at initial year

K	   r	  sigma	rho	  p0
400	 0.5	0.1  	1	    0.7
0.5	 0.5	0.5	  0.001	0.001

The second row relates to the coefficient of variation of each parameter estimation around the expected value (2nd row). A very low value indicates that the parameter is fixed (e.g. 0.01)
