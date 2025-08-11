Rafikisol

What is it?

Rafikisol is a package with helper functions mainly for GEE within R wuch as boostrapping models and categorical downscaling as well as categorical depth functions. We are also
working on many other project. However, many of the functions are in the beta version and if problems arise please report them. The datasets go with a paper that is under review
and hope to come out with more datasets after recieving permission. 

#Example: Categorical splines

#install
remotes::install_github("rafikisol/rafikisol")

#Get data
data(Neno)

#Look at data
head(Neno)

#Run spline
spl = catSpline(Neno, class.var = "Class", type = "class", beta = 1, lam = 0.01)

#See harmonised texture classes
spl$harmonised

#See associated uncertainties
spl$uncertainties

#Evaluate spline (Kappa, accuracy, confusion index, etc.)
evalSpline(spl)

#END
