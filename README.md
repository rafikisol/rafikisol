Rafikisol

What is it?

Rafikisol is a package with helper functions mainly for GEE within R wuch as boostrapping models and categorical downscaling as well as categorical depth functions. We are also
working on many other project. However, many of the functions are in the beta version and if problems arise please report them. The datasets go with a paper that is under review
and hope to come out with more datasets after recieving permission. 

Example: Categorical splines

install

remotes::install_github("rafikisol/rafikisol")

data(Neno)

head(Neno)

spl = holSpline(Neno, class.var = "Class", type = "class", lam = 0.01)

spl$harmonised

spl$uncertainties

evalSpline(spl)

#END
