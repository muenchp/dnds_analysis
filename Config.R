
#### path to precalculated files ####

# #N and #S data
dnds_folder <- "example_data/dnds_values/"
dnds_file_ending <- ".samplePooled.dnds"

# JS data
js_folder <- "example_data/js_score/"
js_file_ending <- ".js"

# Gap data
gap_folder <- "example_data/gap_percent//"
gap_file_ending <- ".js"

# Membrane prediction data
membrane_folder <- "example_data/membrane_prediction/"
membrane_file_ending <- ".tm"

# repeat informations
repeat_folder <- "example_data/repeat_information/"
repeat_file_ending <- ".repeats"


#### values ####
addsmallnumber <- 0.5 # small number to prevent division by zero
window_size <- 10 # window size for sliding window
gapthres <- 0.6 # threshold for gap percentage