# R script that runs AnaCoDa on an input genome.
# Usage: Rscript run_anacoda accession
# where 'accession' is an NCBI accession for which a genome file is available in the current directory

library(AnaCoDa)
library(VGAM)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE) # get command-line input arguments

accession <- args[1] # get accession number

genome <- initializeGenomeObject(file = paste(accession, ".fasta", sep = ''))
parameter <- initializeParameterObject(genome = genome, sphi = 1, num.mixtures = 1, gene.assignment = rep(1, length(genome))) 
mcmc <- initializeMCMCObject(samples = 5000, thinning = 10, adaptive.width=50) 
model <- initializeModelObject(parameter = parameter, model = "ROC") 
runMCMC(mcmc = mcmc, genome = genome, model = model) # this takes a loooong time (like 20 minutes per genome)

# save intermediate files
writeParameterObject(parameter = parameter, file = paste(accession, "-parameter_out.Rda", sep = ''))
writeMCMCObject(mcmc = mcmc, paste(accession, "-mcmc_out.Rda", sep = ''))

# get expression and save
estimatedExpression <- getExpressionEstimates(parameter, 1:length(genome), 1000)
expression <- data.frame(estimatedExpression)
expression$names <- (getNames(genome))
write_csv(expression, paste(accession, "-expression.csv", sep = ''))
