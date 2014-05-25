library(zoo)
library(intervals)

source("Config.R") 
source("ProcessData.R") # collection of data processing functions
source("createPlots.R") # collection of plotting functions


# get total #N and #D value from whole sample
global_values <- globalDnDs(dnds_folder,gap_folder)

createSequencePlot("TIGR00001", global_values$N, global_values$S, global_values$dnds_window) 