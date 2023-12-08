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
# If you have had problems to install the `rstudioapi` package,
# you can find the path to the directory where this script by right clicking
# the name of this script file and selecting "Copy Path". Replace <path> below
# with the path that you have copied (remember to remove the name of the 
# R script from the path, but keep the last "/"!), uncomment the line below,
# and run the command!
#
# setwd('<path>')
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
# X will be an integer (i.e., 9). Subtract two to this number (i.e., 9-2=7)
# and this will be your number of divergence times that are parameters of the
# MCMC. Please modify the number below so it fits to the dataset you are using. 
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
delcol_prior <- 1 # There is one column: mu

# Path to the directory where the subdirectories where each chain ran for each
# dataset are kept. In this case, we will have the path to the directory where
# the analyses when sampling from the prior took place (i.e., this is the path
# to the `CLK` directory that contains the subdirectories from `1` to `5` as
# we ran 5 chains in this example).
## PRIOR ##
path_prior_fosscal <- paste( home_dir, "sum_analyses/prior/fosscal/", sep = "" )
path_prior_seccal  <- paste( home_dir, "sum_analyses/prior/seccal/", sep = "" )
all_paths          <- c( path_prior_fosscal, path_prior_seccal )

#--------------#
# ANALYSE DATA #
#--------------#

# Define object names
num_dirs      <- num_chains
delcol        <- delcol_prior
path          <- all_paths
num_divt      <- num_divt
node_calib    <- calib_nodes
dataset       <- c( "fosscal_CLK", "seccal_CLK" )
perc          <- perc
def_samples   <- def_samples
prior         <- TRUE
out_dat       <- all_paths
time_unit     <- 100
out_file_pdf  <- c( "Convergence_plot_prior_fosscal",
                    "Convergence_plot_prior_seccal" )
out_title_pdf <- c( "CLK | fosscal", "CLK | seccal" )
th            <- 0.1
# Generate a list to keep all the returned objects
sum_CLK_QC          <- vector( "list", length( all_paths ) )
names( sum_CLK_QC ) <- dataset
# Now, run the `QC_conv` function in a loop
for( i in 1:length( dataset ) ){
  # NOTE: Only vectors with more than one element rely on `i`
  sum_CLK_QC[[ i ]] <- QC_conv( num_dirs = num_dirs, delcol = delcol, 
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
save( file = paste( home_dir, "out_RData/sum_CLK_QC.Rdata", sep = "" ),
      sum_CLK_QC )

# Print out sum stats in the script as a log

## [[ FOSSCAL - CLK ]] ####
sum_CLK_QC$fosscal_CLK$ESS_results
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
# Med.    59564    59737
# Min.    58964    58909
# Max.    59818    60180
# 
# $minRhat
# [1] 1.000003
# 
# $maxRhat
# [1] 1.000091
# 
# $total_samples
# [1] 120006      7
length( sum_CLK_QC$fosscal_CLK$not_conv_nodes ) # 0
sum_CLK_QC$fosscal_CLK$num_chains # 6

## [[ SECCAL - CLK ]] ####
sum_CLK_QC$seccal_CLK$ESS_results
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
# Med.    59714    58657
# Min.    57716    56455
# Max.    59982    60035
# 
# $minRhat
# [1] 0.9999685
# 
# $maxRhat
# [1] 1.000117
# 
# $total_samples
# [1] 120006      7
length( sum_CLK_QC$seccal_CLK$not_conv_nodes ) # 0
sum_CLK_QC$seccal_CLK$num_chains # 6

