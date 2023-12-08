#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------------------------------------#
# LOAD PACKAGES, FUNCTIONS, AND SET ENVIRONMENT #
#-----------------------------------------------#
# Load package needed to automatically find the path to the "scripts"
# directory
library( rstudioapi )
scripts_dir   <- gsub( pattern = "scripts..*", replacement = "scripts/",
                       x = getActiveDocumentContext()$path )
setwd( scripts_dir )
# Load the file with all the functions used throughout this script
source( file = "../../../src/Functions.R" )
# Run in-house function to set home directory and output directory for ESS
# and convergence tests
# NOTE: It will create a directory called `plots` and another called
# `ESS_and_chains_convergence` inside the `analyses` directory if you have
# not created them yet
home_dir      <- set_homedir()$home
outchecks_dir <- set_homedir()$ESS
# By now, set the working directory to `home_dir`
setwd( home_dir )

#--------------------------------#
# DEFINE USER'S GLOBAL VARIABLES #
#--------------------------------#
# First, we will define the global variables that match the settings in our 
# analysis.

# 1. Number of chains
num_chains <- 6

# 2. Number of divergence times that have been estimated. One trick to find
# this out quickly is to subtract 1 to the number of species. In this case,
# there are 8 taxa, so the number of internal nodes is `n_taxa-1=8-1=7`.
# Another way to verify this is by opening the `mcmc.txt` file and check the
# header. The first element after `Gen` will have the format of `t_nX`, where
# X will be an integer (i.e., 9). Subtract two to this number and this will
# be your number of divergence times that are parameters of the MCMC. Please
# modify the number below so it fits to the dataset you are using. 
num_divt <- 7

# 3. Number of samples that you specified in the `MCMCtree` control file to 
# collect. NOTE: you may have not collect them all, but do not worry!
def_samples <- 20000

# 4. Quantile percentage that you want to set By default, the variable below is 
# set to 0.975 so the 97.5% and 2.5% quantiles. If you want to change this,
# however, just modify the value.
perc <- 0.975

# 5. Load a semicolon-separated file with info about calibrated nodes. Note that
# each column needs to be separated with semicolons and an extra blank line
# after the last row with calibration information needs to be added (i.e., files
# need to have an extra blank line so R does not complain when reading them). 
# If you add a header, please make sure you name the column elements as 
# `Calib;node;Prior`. If not, the R function below will deal with the header. 
# An example of the format you need to follow to summarise the calibration info
# for each node is the following:
#
# ```
# Calib;node;Prior
# ex_n5;5;ST(5.8300,0.0590,0.1120,109.1240)
# ex_n7;7;B(4.1200,4.5200,0.0250,0.0250)
#
# ```
#
# The first column should have the name of the calibration (e.g., Afrotheria, 
# Laurasiatheria, etc.) as it will help you identify which plot belongs to which
# calibration. The second column is the node used in MCMCtree. The third column
# is the calibration used for that node in MCMCtree format.
# 
# [[ NOTES ABOUT ALLOWED CALIBRATION FORMATS]]
#
# Soft-bound calibrations: 
#  E.g.1: A calibration with a minimum of 0.6 and a maximum of 0.8 would with the 
#         default tail probabilities would have the following equivalent formats:
#         >> B(0.6,0.8) | B(0.6,0.8,0.025,0.025)
#  E.g.2: A calibration with a minimum of 0.6 and a maximum of 0.8 would with the 
#         pL=0.001 and pU=0.025 would have the following format. Note that, whenever
#         you want to modify either pL or pU, you need to write down the four 
#         parameters in the format of "B(min,max,pL,pU)":
#         >> B(0.6,0.8,0.001,0.025)
#
# Lower-bound calibrations: 
#  E.g.1: A calibration with a minimum of 0.6 and the default parameters for
#         p = 0.1, c = 1, pL = 0.025:
#         >> L(0.6) | L(0.6,0.1,1,0.025)
#  E.g.2: A calibration with a hard minimum at 0.6, and so pL = 1e-300. Note that,
#         whenever you want to modify either pL or pU, you need to write down the four 
#         parameters in the format of "L(min,p,c,pL)":
#         >> L(0.6,0.1,1,1e-300)
#
# Upper-bound calibrations: 
#  E.g.1: A calibration with a maximum of 0.8 and the default parameters for pU = 0.025:
#         >> U(0.8) | U(0.8,0.025)
#  E.g.2: A calibration with a hard maximum at 0.8, and so pU = 1e-300. Note that,
#         if you want to modify pU, you need to write down the two
#         parameters in the format of "U(max,pU)":
#         >> U(0.8,1e-300)
#
# ST distributions: 
#  The format accepted has four parameters: xi (location, mean root age), 
#  omega (scale), alpha (shape), nu (df). Accepted format: 
#  >> ST(5.8300,0.0590,0.1120,109.1240)
#
# SN distributions: 
#  The format accepted has three parameters: xi (location, mean root age), 
#  omega (scale), alpha (shape). Accepted format: 
#  >> SN(5.8300,0.0590,0.1120)  
#
#
# The next command executes the `read_calib_f` in-house function, which reads
# your input file (semicolon-separated files). The path to the directory where
# you have saved this file is te argument `main_dir` needs. The argument
# `f_names` requires the file name. Argument `dat` requires a character vector
# to label your dataset. If your input file has a header, please keep
# `head_avail = TRUE`. Otherwise, change this to FALSE.
calib_nodes <- read_calib_f( main_dir = c( paste( home_dir, "calib_files/",
                                                  sep = "" )),
                             f_names = c( "Calibnodes_cal1.csv",
                                          "Calibnodes_cal2.csv" ),
                             dat = c( "bears_fosscal", "bears_seccal" ),
                             head_avail = TRUE )

# 5. Number of columns in the `mcmc.txt` that are to be deleted as they do not 
# correspond to sample values for divergence times (i.e., the entries are not 
# names following the format `t_nX`). To figure out this number quickly, you 
# can open the `mcmc.txt` file, read the header, and count the number of `mu*`
# and `sigma2*` elements. Do not count the `lnL` value when looking at `mcmc.txt`
# files generated when sampling from the posterior -- this is automatically
# accounted for in the in-house R functions that you will subsequently use. 
# E.g., assuming an MCMC ran under a relaxed-clock model with no partitions, 
# we would see `mu` and `sigma2` columns. Therefore, the variable would be set 
# to `delcol = 2`. Please modify the values below according to your dataset for 
# the `mcmc.txt` file generated when sampling from the prior (`delcol_prior`) 
# and when sampling from the posterior (`delcol_posterior`).
delcol_post   <- 2 # There are two column: mu, sigma2

# Path to the directory of each data alignment.
# In this case, we will have the path to the directory where
# the analyses when sampling from the prior took place (i.e., this is the path
# to the `CLK` directory that contains the subdirectories from `1` to `5` as
# we ran 5 chains in this example).
## POSTERIOR ##
path_GBM_conc_fosscal <- paste( home_dir,"sum_analyses/posterior/conc/fosscal/GBM/",
                          sep = "" )
path_GBM_conc_seccal <- paste( home_dir,"sum_analyses/posterior/conc/seccal/GBM/",
                          sep = "" )
path_GBM_part_fosscal <- paste( home_dir,"sum_analyses/posterior/part/fosscal/GBM/",
                               sep = "" )
path_GBM_part_seccal <- paste( home_dir,"sum_analyses/posterior/part/seccal/GBM/",
                                sep = "" )
path_ILN_conc_fosscal <- paste( home_dir,"sum_analyses/posterior/conc/fosscal/ILN/",
                                sep = "" )
path_ILN_conc_seccal <- paste( home_dir,"sum_analyses/posterior/conc/seccal/ILN/",
                               sep = "" )
path_ILN_part_fosscal <- paste( home_dir,"sum_analyses/posterior/part/fosscal/ILN/",
                                sep = "" )
path_ILN_part_seccal <- paste( home_dir,"sum_analyses/posterior/part/seccal/ILN/",
                               sep = "" )

## Add only those paths in the following variable that your team has analysed!
## Remove from the command below those that do not correspond to the analyses 
## you have run!
all_paths <- c( path_GBM_conc_fosscal, path_GBM_conc_seccal,
                path_GBM_part_fosscal, path_GBM_part_seccal,
                path_ILN_conc_fosscal, path_ILN_conc_seccal,
                path_ILN_part_fosscal, path_ILN_part_seccal )

#--------------#
# ANALYSE DATA #
#--------------#

# Define object names
num_dirs      <- num_chains
delcol        <- delcol_post
path          <- all_paths
num_divt      <- num_divt
## Remove from the options below those referring to the datasets
## that you have not analysed!
dataset       <- c( "GBM_conc_fosscal", "GBM_conc_seccal",
                    "GBM_part_fosscal", "GBM_part_seccal",
                    "ILN_conc_fosscal", "ILN_conc_seccal",
                    "ILN_part_fosscal", "ILN_part_seccal" )
## Remove from the options below those referring to the datasets
## that you have not analysed!
node_calib    <- list( calib_nodes$bears_fosscal, calib_nodes$bears_seccal,
                       calib_nodes$bears_fosscal, calib_nodes$bears_seccal,
                       calib_nodes$bears_fosscal, calib_nodes$bears_seccal,
                       calib_nodes$bears_fosscal, calib_nodes$bears_seccal ) 
names( node_calib ) <- dataset
perc          <- perc
def_samples   <- def_samples
prior         <- FALSE
out_dat       <- all_paths
time_unit     <- 100
## Remove from the options below those referring to the datasets
## that you have not analysed!
out_file_pdf  <- c( "Convergence_plot_post_GBM_conc_fosscal",
                    "Convergence_plot_post_GBM_conc_seccal",
                    "Convergence_plot_post_GBM_part_fosscal",
                    "Convergence_plot_post_GBM_part_seccal",
                    "Convergence_plot_post_ILN_conc_fosscal",
                    "Convergence_plot_post_ILN_conc_seccal",
                    "Convergence_plot_post_ILN_part_fosscal",
                    "Convergence_plot_post_ILN_part_seccal" )
## Remove from the options below those referring to the datasets
## that you have not analysed!
out_title_pdf <- c( "GBM and all CPs | fosscal",
                    "GBM and all CPs | seccal",
                    "GBM and 12CPs | fosscal",
                    "GBM and 12CPs | seccal",
                    "ILN and all CPs | fosscal",
                    "ILN and all CPs | seccal",
                    "ILN and 12CPs | fosscal",
                    "ILN and 12CPs | seccal" )
th            <- 0.1
# Generate a list to keep all the returned objects
sum_post_QC          <- vector( "list", length( all_paths ) )
names( sum_post_QC ) <- dataset
# Now, run the `QC_conv` function in a loop
for( i in 1:length( dataset ) ){
  # NOTE: Only vectors with more than one element rely on `i`
  sum_post_QC[[ i ]] <- QC_conv( num_dirs = num_dirs, delcol = delcol, 
                                 path = path[i], num_divt = num_divt, 
                                 node_calib = node_calib[[ i ]],
                                 dataset = dataset[i],
                                 perc = perc, def_samples = def_samples,
                                 prior = prior, out_dat = out_dat[i],
                                 time_unit = time_unit, 
                                 out_file_pdf = out_file_pdf[i],
                                 out_title_pdf = out_title_pdf[i], 
                                 th = th )
}

# Save object!
if( ! dir.exists( paste( home_dir, "out_RData", sep = "" ) ) ){
  dir.create( paste( home_dir, "out_RData", sep = "" ) ) 
}
save( file = paste( home_dir, "out_RData/sum_post_QC", sep = "" ),
      sum_post_QC )


# Print out sum stats in the script as a log

## [[ FOSSCAL - CONC - GBM ]] ####
sum_post_QC$GBM_conc_fosscal$ESS_results
# $median
# [1] 20001
# 
# $min
# [1] 20001
# 
# $max
# [1] 20001
# 
# $num_samples_for_ESS
# [1] 120006
# 
# $ESS
# Tail-ESS Bulk-ESS
# Med.    59650    60334
# Min.    59074    57192
# Max.    60448    60832
# 
# $minRhat
# [1] 0.999959
# 
# $maxRhat
# [1] 1.000006
# 
# $total_samples
# [1] 120006      7
length( sum_post_QC$GBM_conc_fosscal$not_conv_nodes ) # 0
sum_post_QC$GBM_conc_fosscal$num_chains # 6


## [[ SECCAL - CONC - GBM ]] ####
sum_post_QC$GBM_conc_seccal$ESS_results
# $median
# [1] 20001
# 
# $min
# [1] 20001
# 
# $max
# [1] 20001
# 
# $num_samples_for_ESS
# [1] 120006
# 
# $ESS
# Tail-ESS Bulk-ESS
# Med.    59272    59820
# Min.    56769    59067
# Max.    60318    60345
# 
# $minRhat
# [1] 1.000011
# 
# $maxRhat
# [1] 1.000035
# 
# $total_samples
# [1] 120006      7
length( sum_post_QC$GBM_conc_seccal$not_conv_nodes ) # 0
sum_post_QC$GBM_conc_seccal$num_chains # 6


## [[ FOSSCAL - PART - GBM ]] ####
sum_post_QC$GBM_part_fosscal$ESS_results
# $median
# [1] 20001
# 
# $min
# [1] 20001
# 
# $max
# [1] 20001
# 
# $num_samples_for_ESS
# [1] 120006
# 
# $ESS
# Tail-ESS Bulk-ESS
# Med.    58916    59470
# Min.    56848    57976
# Max.    59609    59782
# 
# $minRhat
# [1] 0.9999801
# 
# $maxRhat
# [1] 1.000022
# 
# $total_samples
# [1] 120006      7
length( sum_post_QC$GBM_part_fosscal$not_conv_nodes ) # 0
sum_post_QC$GBM_part_fosscal$num_chains # 6


## [[ SECCAL - PART - GBM ]] ####
sum_post_QC$GBM_part_seccal$ESS_results
# $median
# [1] 20001
# 
# $min
# [1] 20001
# 
# $max
# [1] 20001
# 
# $num_samples_for_ESS
# [1] 120006
# 
# $ESS
# Tail-ESS Bulk-ESS
# Med.    59652    59931
# Min.    59118    57787
# Max.    60120    60764
# 
# $minRhat
# [1] 0.999986
# 
# $maxRhat
# [1] 1.000168
# 
# $total_samples
# [1] 120006      7
length( sum_post_QC$GBM_part_seccal$not_conv_nodes ) # 0
sum_post_QC$GBM_part_seccal$num_chains # 6


## [[ FOSSCAL - CONC - ILN ]] ####
sum_post_QC$ILN_conc_fosscal$ESS_results
# $median
# [1] 20001
# 
# $min
# [1] 20001
# 
# $max
# [1] 20001
# 
# $num_samples_for_ESS
# [1] 120006
# 
# $ESS
# Tail-ESS Bulk-ESS
# Med.    58804    57466
# Min.    58099    56795
# Max.    59886    60190
# 
# $minRhat
# [1] 1.000006
# 
# $maxRhat
# [1] 1.000102
# 
# $total_samples
# [1] 120006      7
length( sum_post_QC$ILN_conc_fosscal$not_conv_nodes ) # 0
sum_post_QC$ILN_conc_fosscal$num_chains # 6

## [[ SECCAL - CONC - ILN ]] ####
sum_post_QC$ILN_conc_seccal$ESS_results
# $median
# [1] 20001
# 
# $min
# [1] 20001
# 
# $max
# [1] 20001
# 
# $num_samples_for_ESS
# [1] 120006
# 
# $ESS
# Tail-ESS Bulk-ESS
# Med.    59230    59405
# Min.    57962    58307
# Max.    60101    59787
# 
# $minRhat
# [1] 0.9999349
# 
# $maxRhat
# [1] 1.000129
# 
# $total_samples
# [1] 120006      7
length( sum_post_QC$ILN_conc_seccal$not_conv_nodes ) # 0
sum_post_QC$ILN_conc_seccal$num_chains # 6


## [[ FOSSCAL - PART - ILN ]] ####
sum_post_QC$ILN_part_fosscal$ESS_results
# $median
# [1] 20001
# 
# $min
# [1] 20001
# 
# $max
# [1] 20001
# 
# $num_samples_for_ESS
# [1] 120006
# 
# $ESS
# Tail-ESS Bulk-ESS
# Med.    59113    59350
# Min.    57932    58930
# Max.    59576    60175
# 
# $minRhat
# [1] 0.9999376
# 
# $maxRhat
# [1] 1.00003
# 
# $total_samples
# [1] 120006      77
length( sum_post_QC$ILN_part_fosscal$not_conv_nodes ) # 0
sum_post_QC$ILN_part_fosscal$num_chains # 6

## [[ SECCAL - PART - ILN ]] ####
sum_post_QC$ILN_part_seccal$ESS_results
# $median
# [1] 20001
# 
# $min
# [1] 20001
# 
# $max
# [1] 20001
# 
# $num_samples_for_ESS
# [1] 120006
# 
# $ESS
# Tail-ESS Bulk-ESS
# Med.    59156    59168
# Min.    58146    56946
# Max.    59703    59919
# 
# $minRhat
# [1] 1.000021
# 
# $maxRhat
# [1] 1.000199
# 
# $total_samples
# [1] 120006      7
length( sum_post_QC$ILN_part_seccal$not_conv_nodes ) # 0
sum_post_QC$ILN_part_seccal$num_chains # 6
