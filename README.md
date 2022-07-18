# Aplomado-Falcon-IPM

[![DOI](https://zenodo.org/badge/291721348.svg)](https://zenodo.org/badge/latestdoi/291721348)

Supplemental materials for: B. W. Rolek, B. Pauli, A. Macias-Duarte, B. Mutch, P. Juergens, T. Anderson, C. Parish, J. Johnson, B. Millsap, C. J. W. McClure. 2022. Long-term Demography of a Reintroduced and Isolated Population of Endangered Falcons. Global Ecology and Conservation.
 
 Contact author: brianrolek at gmail.com.
 
Metadata associated with the file "data\data-7states.Rdata". Once loaded into R, this file contains the object "datl", a list containing data for input into the JAGS model. Data fields include: 

+ y = observed states: 1 = seen first year (age 0), 2 = seen nonbreeder, 3 = seen breeder, 4 = recovered dead, 5 = not seen. Data are a matrix of dimensions nind x n.yr.
+ first = the first year of capture and banding.A vector of length nind.
+	n.yr = number of years. A scalar.
+	nind = number of individuals. A scalar.
+	hacked = index indicating whether an individual's origin. 1 = wild-hatched, 2 = hacked. A vector of length nind.		
+ sex = index indicating whether an individual was a male or female. 
		        1 = female, 2 = male. A vector of length nind.
+ l.countBM = log(counts) of adult breeder males each year. Zeroes evaluated to "Inf" so they are imputed with a small constant of -6. A vector of length n.yr.
+ l.countFM = log(counts) of adult nonbreeder males (floaters) each year. See l.countBM.
+	l.countJM = log(counts) of first-year males (juveniles) each year. See l.countBM.
+	countBM = counts of adult breeder males each year. A vector of length n.yr.
+	countFM = counts of adult nonbreeder males each year. A vector of length n.yr.
+	countJM = counts of first-year males each year. A vector of length n.yr.
+	countOM = counts of wild-hatched first-year males each year. A vector of length n.yr.
+	countBF = counts of adult breeder females each year. A vector of length n.yr.
+	countFF = counts of adult nonbreeder females each year. A vector of length n.yr.
+	countJF = counts of first-year females each year. A vector of length n.yr.
+	countOF = counts of wild-hatched first-year females each year. A vector of length n.yr.
+	aug = the number of first-year captive-bred birds translocated for each year.  A vector of length n.yr.
+	manage = an index for each year indicating whether improvements to nest boxes or habitat had occurred. 1 = pre-management <2012, 2= post-management >=2012.  A vector of length n.yr. 
+	prod = productivity for each spatial territory and breeding pair. A vector of length "number of nests monitored".
+	pp = total productivity for each year. Referred to as "J" in manuscript. A vector of length n.yr.
+	year.p = an index of years for prod. A vector of length n.ter.
+	ter = an index of spatial territories for prod. A vector of length n.ter.
+	K = the number of observed nests. Scalar. 
+	n.ter = the total number of territories.
+	effort = an index for each year indicating whether survey effort was 1 = low, or 2 = high. A vector of length n.yr.
+	z =  data of known states. calculated using the function known.state.ms() in the model file. 1 = alive first year, 2 = alive nonbreeder, 3 = alive breeder, 4 = dead, 5 = dead not recovered/long dead, 6 = emigrated alive, 7 = emigrated dead. States 5, 6, and 7 are unobserved and only occur in the CJS model that includes emigration. Data are a matrix of dimensions nind x n.yr.
