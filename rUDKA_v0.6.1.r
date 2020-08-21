# Type:	        R-Program

# Title:        R framework for a Unified Dispersal Kernel Analysis - rUDKA

# Aims:         Function Fitting of Dispersal Kernels

# Version:      0.6.1

# Author/s:     Dr. Reinhard Klenke & Dr. Rebecca Thier-Lange

# Institution:  Helmholtz-Centre for Environmental Research
# Department:   Department of Conservation Biology
# Street:       Permoser Str. 15
# Town:         04318 Leipzig
# Country:      G E R M A N Y
# URL:          http://www.ufz.de


# Email:        reinhard.klenke [at] ufz.de
#               rebecca.lange [at] ufz.de

# Abstract:
# This program is a compilation of methods used to fit different functions
# to organism dispersal kernels

# Citation:
# Klenke A.R., Thier-Lange R (2020): rUDKA - An R framework for a Unified Dispersal
# Kernel Analysis. R program code. Available at https://github.com/popecologist/rUDKA. 

# Keywords:    Animal Dispersal, Dispersal Kernel, Growth Functions

# We used te grofit library for fitting:
# Kahm M, Hasenbrink G, Lichtenberg-Frate H, Ludwig J, Kschischo M (2010): grofit: Fitting
# Biological Growth Curves with R. Journal of Statistical Software, 33(7), 1-21. 
# URL http://www.jstatsoft.org/v33/i07/."

# License:
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# Please inform us (see email and address above) if you find bugs 
# and mistakes.



# History:
# -------------------------------------------------------------------- 
# Changing title from "Function Fitting of Dispersal Kernels" to 
# "R framework for a Unified Dispersal Kernel Analysis"
# Adding "Aims: Function Fitting of Dispersal Kernels"
# Adding recommendation for citation
# Formatting some comments
# Adding recommendations for the content of the header of the output
# 10.08.2020 13:27:00
# bootstrapping was commented out because it caused some problems 
# with spline fit
# 27.11.2012 19:17:06
# Some small improvement regarding legend positioning 
# 20.11.2012 00:11:48
# Some small improvements regarding data handling (e.g. na.omit())
# 19.11.2012 23:41:21
# First running version of the program. Only growth functions are used.
# 19.11.2012 22:28:59
# -------------------------------------------------------------------- #

# Clean the workspace from all objects
rm (list = ls ())

# Call the garbage collector and free memory
gc ()

# load libraries
library (nlstools)
library (grofit)
library (growthrates)
library (car)



# -------------------------------------------------------------------- #
# Declaration of special functions used in the analysis

# useful for placement of legends
goldenCut <- function (x, pos = c ("l","r")) {
	ifelse (pos == "l", goldenCut <- max (x) * 1/exp (1), goldenCut <- max (x) * (1- 1/exp (1)))
	return (goldenCut)
}


# "saveGraphs"
# Saves graphics in both formats WMF (pdf for unix/MacOS) and PNG in two
# separate directories (outWMF and outHTML) which have to be declared
# in the main program!
 
saveGraphs <- function (Var)
     {
 		if(.Platform$OS.type == "windows") {
 			savePlot (filename = file.path (outWMF, Var),
               	type="wmf",
               	device=dev.cur())
            savePlot (filename =file.path (outHTML, Var),
               type="png",
               device=dev.cur())
 					} else {
 			dev.copy2pdf (file=file.path (outWMF, paste(Var,"pdf",sep=".")))
 					}      
     dev.off ()
     }

# End of declaration of special functions used in the analysis
# -------------------------------------------------------------------- #



# Declaration of fitting functions

# Exponential distribution 1
class (f_exp1 <- nrInd ~ a * exp (-b * distances))

# Exponential distribution 2
class (f_exp2 <- nrInd ~ exp (a - b * distances^2))

# Sigmoidal distribution 1  from Guy used in FunConn
class (f_sig1 <- nrInd ~ (exp (a * (distances - b)) - exp(-a * (distances - b))) / (exp(a * (distances - b)) + exp(-a * (distances - b))))

# Sigmoidal distribution 2
class (f_sig2 <- nrInd ~ 1 - exp (-a * exp(-b * distances)))

# Inverse power function
class (f_ipf1 <- nrInd ~ a * distances^(-b))

# Inverse power function
class (f_ipf2 <- nrInd ~ a / (1+distances^(-b)))

# Fat tailed dispersal kernel
class (f_fatT <- nrInd ~ 1 / (1 + a * distances^-b))

#            y~ ((A^x / (factorial(x)) * exp(-A)) * sum(results[,2])
#            y~ ((A^x / (sqrt(2 * Pi * x) * (x / exp)^x)) * exp(-A)) * sum(results[,2])     substituting 'factorial(x)'

# Poisson distribution
class (f_pois <- nrInd ~ dpois (distances, lambda))

# Gamma distribution
class (f_gamma <- nrInd ~ dgamma (distances, shape, scale))

# Polynom
class (f_poly <- nrInd ~ a + a * distances + b * distances^2 + c * distances^3)


# Poisson distribution
# The Poisson distribution arises in connection with Poisson processes. 
# It applies to various phenomena of discrete properties (that is, 
# those that may happen 0, 1, 2, 3, ... times during a given period 
# of time or in a given area) whenever the probability of the phenomenon 
# happening is constant in time or space. Examples of events that may be 
# modelled as a Poisson distribution include:
#
# ...
#    The number of mutations in a given stretch of DNA after a certain amount of radiation.
#    The proportion of cells that will be infected at a given multiplicity of infection.
#
# lambda is approximately the median of the distribution
#
# scale <- number of counts at the median of the distribution / 
# (sum (counts) - number of counts at the median of the distribution)
#
# scale <- number of counts at the median of the distribution / scale

class (f_pois <- NrInd ~ (lambda^distances / factorial (distances) * exp (-lambda)))

# Gamma distribution
class (f_dgamma <- NrInd~dgamma (distances, shape = p, scale = b) * y)



# -------------------------------------------------------------------- #

# preset of parameters
a <- 0.01
b <- 0.01
c <- 0.05
d <- 0.05

lambda <- 10
rate = b
scale = 1/rate

# -------------------------------------------------------------------- #


# Define report titles
ReportTitle <- "Dispersal kernel Analysis"

# Define working directory and filenames

# WorkDir <- "drive:/directory"
WorkDir <- getwd ()
# WorkDir <- "rUDKA-Dispersal"

OUTFile <- "results.out"      # Output of the program

# Define subdirectories for output
OutTXT <- "OutTXT"
OutHTM <- "OutHTM"
OutWMF <- "OutWMF"
InCSV <- "InCSV"

# Set working directory
setwd (WorkDir)

# Create output directories
dir.create (file.path (WorkDir,OutTXT))
dir.create (file.path (WorkDir,OutWMF))
dir.create (file.path (WorkDir,OutHTM))

# Declare variables
outTXT <- file.path (WorkDir,OutTXT)
outWMF <- file.path (WorkDir,OutWMF)
outHTML <- file.path (WorkDir,OutHTM)
inCSV <-  file.path (WorkDir,InCSV)

# Define output files
outFileTXT <- paste (outTXT,
                     OUTFile,
                     sep="/")

# Split output to the console and a file
outResultsTXT <- file (outFileTXT,
                       "w+")

sink (outResultsTXT,
      split= TRUE,
      append = FALSE)


# read file list
inFiles <- dir (path = inCSV, pattern = ".csv")

# calculate file list length
k <- length (inFiles)

# define plot title based on species and analysis typ = input file name
plotTitle <- vector (mode = "character", length = k)

# initialize variables
resultTableRow <- 0
resultTable <- matrix (data = NA, nrow = k, ncol = 37)

# preparing indices to get access on data in objects
RealData <- 1
NormData <- 2

# Use the following FOR loop OR set i and k to the index you want for stepwise calculations
# skip the FOR declaration and the closing bracket at the end of the loop for testing
for (i in 1:k) {

# set i = 1 for testing and remove the "#" otherwise comment the following line out!
# i = 1

    dispersal.kernel <- read.table(paste(inCSV, inFiles [i], sep="/"), 
                                   sep = ";", dec = ",", header = TRUE)
    dispersal.kernel = na.omit(dispersal.kernel)
    names (dispersal.kernel) <- c("distances", "nrInd")
    dispersal.kernel$cumNrInd <- cumsum(dispersal.kernel$nrInd)
    dispersal.kernel$stdNrInd <- dispersal.kernel$nrInd / max (dispersal.kernel$cumNrInd)
    dispersal.kernel$stdCumNrInd <- dispersal.kernel$cumNrInd / max (dispersal.kernel$cumNrInd)
    dispersal.kernel$stdDist <- dispersal.kernel$distances / max (dispersal.kernel$distances)
  
    opar <- par (mfrow = c(2,1), cex.main = 1)
    
    plotTitle [i] <- substring (inFiles [i], 1, (nchar (inFiles [i])-4))
    plot (dispersal.kernel$nrInd~dispersal.kernel$distances,
          main = plotTitle [i], xlab = "Distance", ylab = "Frequency",
          type = "h", lwd = 10, col = "darkgrey", lend = 2)

    data.smoothed <- smooth.spline (dispersal.kernel$distances, dispersal.kernel$nrInd)
    
    lines (data.smoothed$y ~ data.smoothed$x, col = "salmon", lwd = 3)

    plot (dispersal.kernel$cumNrInd~dispersal.kernel$distances,
          main = plotTitle [i], xlab = "Distance", ylab = "Cumulative Frequency")

    data.smoothed <- smooth.spline (dispersal.kernel$distances, dispersal.kernel$cumNrInd)
    
    lines (data.smoothed$y ~ data.smoothed$x, col = "salmon", lwd = 3)

# insert saveFiles!!!

# Save graph

    CurrentFileName <- paste (plotTitle [i],
                              "descriptive",
                              "freq",
                              sep = "_")
    
    saveGraphs (CurrentFileName)
    
    par (opar)
    par (mfrow = c(2,1), cex.main = 1)

    plot (dispersal.kernel$stdNrInd~dispersal.kernel$stdDist,
          main = plotTitle [i], xlab = "Std.Distance", ylab = "Density",
          type = "h", lwd = 10, col = "darkgrey", lend = 2)

    data.smoothed <- smooth.spline (dispersal.kernel$stdDist, dispersal.kernel$stdNrInd, 
                                  spar=0.5)

    lines (data.smoothed$y ~ data.smoothed$x, col = "salmon", lwd = 3)

    plot (dispersal.kernel$stdCumNrInd~dispersal.kernel$stdDist, 
        main = plotTitle [i], xlab = "Distance", ylab = "Probability")

    data.smoothed <- smooth.spline (dispersal.kernel$stdDist, dispersal.kernel$stdCumNrInd, 
                                    spar= 0.5)
    
    lines (data.smoothed$y ~ data.smoothed$x, col = "salmon", lwd = 3)

# insert saveFiles!!!

# Save graph
    CurrentFileName <- paste (plotTitle [i],
                              "descriptive",
                              "prob",
                              sep = "_")
    
    saveGraphs (CurrentFileName)


# prepare data structure
    dispersal.desc <- rbind (c(plotTitle [i], "Real Data"),
                             c(plotTitle [i], "Normalized Data"))

    dispersal.dose <- rbind (0,0)

    dispersal.data <- rbind (dispersal.kernel$cumNrInd,
                             dispersal.kernel$stdCumNrInd)

    dispersal.data <- cbind (data.frame(dispersal.desc), 
                             data.frame (dispersal.dose),
                             data.frame(dispersal.data))

    dispersal.dist <- data.frame(rbind (dispersal.kernel$distances,
                                        dispersal.kernel$stdDist))

    try (resultsGroFit <- grofit (time = dispersal.dist, 
                                  data = dispersal.data,
                                  ec50 = FALSE,
                                  control = grofit.control(fit.opt = "b",
                                                           # bootstrapping not necessary
                                                           nboot.gc = 0, 
                                                           interactive = FALSE)),
         silent = FALSE
        )

    try (resultsGcFit <- gcFit (time = dispersal.dist, 
                                data = dispersal.data,
                                control=grofit.control(fit.opt="b", 
                                                       interactive = FALSE)),
         silent = FALSE)

    try (summary(resultsGroFit))

    try (parameters <- summary (resultsGroFit$gcFit), 
                                silent = FALSE)


    modelOutput <- resultsGcFit$gcFittedModels[[RealData]]

    resultsNls <- summary (modelOutput$nls)
    residualsNls <- resultsNls$residuals

    xanchor <- goldenCut (modelOutput$raw.time, pos = "l")
    yanchor <- goldenCut (modelOutput$raw.data, pos = "r")

    par (mfrow = c(2,2))

    plot (modelOutput$raw.time, 
          modelOutput$raw.data, 
          xlab = "Distance",
          ylab = "Cumulative Frequency", 
          main = plotTitle [i],
          cex.main = 0.6,
          ylim = c(0, max (modelOutput$raw.data)))

    lines (modelOutput$fit.time, modelOutput$fit.data, lwd = 2, col = "blue")


    bla <- modelOutput$fit.time * parameters$mu.model [RealData]
    bla <- bla + (-parameters$mu.model [RealData] * parameters$lambda.model [RealData])
    
    lines(modelOutput$fit.time, bla, lw = 2, lty = 2, col = "red")
    legend (xanchor, yanchor, legend = c(expression (lambda), c(expression (mu)), "A"), bt = "n")
    legend (xanchor+max(xanchor)*0.2, yanchor,
            legend = c(paste ("=",round (parameters$lambda.model [RealData],3)),
                       paste ("=",round (parameters$mu.model [RealData],3)),
                       paste ("=",round (parameters$A.model [RealData],3))), bt="n")


    modelOutput <- resultsGcFit$gcFittedModels[[NormData]]

    xanchor <- goldenCut (modelOutput$raw.time, pos = "l")
    yanchor <- goldenCut (modelOutput$raw.data, pos = "r")

    plot (modelOutput$raw.time, 
          modelOutput$raw.data, 
          xlab = "Distance",
          ylab = "Probability", 
          main = "Fitted Model",
          ylim = c(0, max (modelOutput$raw.data)))

    lines (modelOutput$fit.time, modelOutput$fit.data, lwd = 2, col = "blue")

    bla <- modelOutput$fit.time * parameters$mu.model [NormData]
    bla <- bla + (-parameters$mu.model [NormData] * parameters$lambda.model [NormData])
    
    lines(modelOutput$fit.time, bla, lw = 2, lty = 2, col = "red")

    legend (xanchor, yanchor, 
            legend = c(expression (lambda), 
                       c(expression (mu)), 
                       "A"), bt="n")
    legend (xanchor+max(xanchor)*0.2, yanchor,
            legend = c(paste ("=",round (parameters$lambda.model [NormData],3)),
                       paste ("=",round (parameters$mu.model [NormData],3)),
                       paste ("=",round (parameters$A.model [NormData],3))), bt="n")

    splineOutput <- resultsGcFit$gcFittedSplines [[RealData]]
    xanchor <- goldenCut (splineOutput$raw.time, pos = "l")
    yanchor <- goldenCut (splineOutput$raw.data, pos = "r")

    plot (splineOutput$raw.time, 
          splineOutput$raw.data, 
          xlab = "Distance",
          ylab = "Cumulative Frequency", 
          main = plotTitle [i],
          cex.main = 0.6,
          ylim = c(0, max (splineOutput$raw.data)))

    lines (splineOutput$fit.time, splineOutput$fit.data, lwd = 2, col = "blue")

    bla <- splineOutput$fit.time * parameters$mu.spline [RealData]
    bla <- bla + (-parameters$mu.spline [RealData] * parameters$lambda.spline [RealData])
    lines(splineOutput$fit.time, bla, lw = 2, lty = 2, col = "red")

    legend (xanchor, yanchor, 
            legend = c(expression (lambda), 
                       c(expression (mu)), 
                       "A"), bt = "n")
    legend (xanchor+max(xanchor)*0.2, yanchor,
            legend = c(paste ("=",round (parameters$lambda.spline [RealData],3)),
                       paste ("=",round (parameters$mu.spline [RealData],3)),
                       paste ("=",round (parameters$A.spline [RealData],3))), bt="n")

    splineOutput <- resultsGcFit$gcFittedSplines [[NormData]]

    xanchor <- goldenCut (splineOutput$raw.time, pos = "l")
    yanchor <- goldenCut (splineOutput$raw.data, pos = "r")

    plot (splineOutput$raw.time, 
          splineOutput$raw.data, 
          xlab = "Distance",
          ylab = "Probability", 
          main = "Fitted Spline",
          ylim = c(0, max (splineOutput$raw.data)))

    lines (splineOutput$fit.time, splineOutput$fit.data, lwd = 2, col = "blue")

    bla <- splineOutput$fit.time * parameters$mu.spline [NormData]
    bla <- bla + (-parameters$mu.spline [NormData] * parameters$lambda.spline [NormData])
    
    lines(splineOutput$fit.time, bla, lw = 2, lty = 2, col = "red")

    legend (xanchor, yanchor, 
            legend = c(expression (lambda), 
                       c(expression (mu)), "A"), bt = "n")
    
    legend (xanchor+max(xanchor)*0.2, yanchor,
            legend = c(paste ("=",round (parameters$lambda.spline [NormData],3)),
                       paste ("=",round (parameters$mu.spline [NormData],3)),
                       paste ("=",round (parameters$A.spline [NormData],3))), bt="n")

  # Save graph
    
    CurrentFileName <- paste (plotTitle [i],
                              "descriptive_fit",
                              sep = "_")

    saveGraphs (CurrentFileName)


# save parameters to a structure
# don't forget to change 1 to i
    
    resultTableRow <- c (plotTitle [i],
                      parameters$used.model [RealData],
                      parameters$A.model [RealData],
                      parameters$stdA.model [RealData],
                      parameters$ci95.A.model.lo [RealData],
                      parameters$ci95.A.model.up [RealData],
                      parameters$mu.model [RealData],
                      parameters$stdmu.model [RealData],
                      parameters$ci95.mu.model.lo [RealData],
                      parameters$ci95.mu.model.up [RealData],
                      parameters$lambda.model [RealData],
                      parameters$stdlambda.model [RealData],
                      parameters$ci95.lambda.model.lo [RealData],
                      parameters$ci95.lambda.model.up [RealData],
                      parameters$integral.model [RealData],
                      parameters$mu.spline [RealData],
                      parameters$lambda.spline [RealData],
                      parameters$A.spline [RealData],
                      parameters$integral.spline [RealData],
                      parameters$used.model [NormData],
                      parameters$A.model [NormData],
                      parameters$stdA.model [NormData],
                      parameters$ci95.A.model.lo [NormData],
                      parameters$ci95.A.model.up [NormData],
                      parameters$mu.model [NormData],
                      parameters$stdmu.model [NormData],
                      parameters$ci95.mu.model.lo [NormData],
                      parameters$ci95.mu.model.up [NormData],
                      parameters$lambda.model [NormData],
                      parameters$stdlambda.model [NormData],
                      parameters$ci95.lambda.model.lo [NormData],
                      parameters$ci95.lambda.model.up [NormData],
                      parameters$integral.model [NormData],
                      parameters$mu.spline [NormData],
                      parameters$lambda.spline [NormData],
                      parameters$A.spline [NormData],
                      parameters$integral.spline [NormData])

resultTable [i,] <- resultTableRow

# uncomment following bracket to make the FOR loop working
}



# naming parameter data.frame

colnames (resultTable) <- c("Species Kernel","Model.used", "A.model","stdA.model","ci95.A.model.lo", "ci95.A.model.up",
                         "mu.model", "stdmu.model", "ci95.mu.model.lo", "ci95.mu.model.up",
                         "lambda.model", "stdlambda.model", "ci95.lambda.model.lo", "ci95.lambda.model.up",
                         "integral.model",
                         "mu.spline", "lambda.spline", "A.spline", "integral.spline",
                         "Model.norm.used",
                         "A.norm.model","stdA.norm.model", "ci95.A.norm.model.lo", "ci95.A.norm.model.up",
                         "mu.norm.model", "stdmu.norm.model", "ci95.mu.norm.model.lo", "ci95.mu.norm.model.up",
                         "lambda.norm.model", "stdlambda.norm.model", "ci95.lambda.norm.model.lo", "ci95.lambda.norm.model.up",
                         "integral.norm.model",
                         "mu.norm.spline", "lambda.norm.spline", "A.norm.spline",
                         "itegral.norm.spline")

# Adding all important information (what, who, which data, when) about your analysis here
header <- c("# Dispersal Kernel Analysis from:",
            timestamp (),
            "# Done by",
            "# Reinhard Klenke",
            "# Data extracted from literature and prepared/revised by",
            "# Rebecca Thier-Lange",
            "# Claudia Guimaraes-Steinicke",
            "# Leonie Bartel",
            "# Jana Bermudez Alvarez", 
            "# Guy Pe'er", 
            "# Reinhard Klenke",
            "# ... more comments ...",
            "#",
            "# Fitting was done with the grofit library described in:",
            "# Kahm M, Hasenbrink G, Lichtenberg-Frate H, Ludwig J, Kschischo M (2010)",
            "# grofit: Fitting Biological Growth Curves with R.",
            "# Journal of Statistical Software, 33(7), 1-21.",
            "# URL http://www.jstatsoft.org/v33/i07/."))


# save data as file
write.table (header,
       file = paste (outTXT,"parameters.csv",sep = "/"),
       quote = FALSE,
       row.names = FALSE,
       col.names = FALSE,
       sep = "\n")

write.table (resultTable,
       file = paste (outTXT,"parameters.csv",sep = "/"),
       sep = ";",
       row.names = FALSE,
       col.names = TRUE,
       append = TRUE)

# save parameters to a file


# End of text file output and closing of files
sink ()

# Clean the workspace from all objects
rm (list = ls ())

# Call the garbage collector to free the computer memory from wasted
# objects

gc ()

# -------------------------------------------------------------------- #

# End of program
