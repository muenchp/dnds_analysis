<<<<<<< HEAD
=======
# autor
>>>>>>> 79c713afd62f09ba7294774892c0778500191683

globalDnDs <- function(path_to_folder_with_dnds_files, path_to_folder_with_gap_files){
  ### get number of synonymous and non-synonymous mutations for a whole sample ###
  # requires: 
  # packages: zoo, 
  # input: path_to_folder_with_dnds_files
  #        path_to_folder_with_gap_files
  # output list object with (1) #N (2) #S (3) #N/#S in sliding window
  
  tigr_list <- list.files(path_to_folder_with_dnds_files) # list of all input files 
  dnds_sample_all <- NULL; i <- 0; dNall <- 0; dSall <- 0 # assign variables
<<<<<<< HEAD
=======
  addsmallnumber <- 0.5 # small number to prevent division by zero
  window_size <- 10 
  gapthres <- 0.6
>>>>>>> 79c713afd62f09ba7294774892c0778500191683
  for (file in tigr_list) {
    # iterate over all families in input folder
    i <- i + 1
    tigr <- substr(tigr_list[i], 1, 9) 
    cat(paste("Processing", tigr, "\n"))
    
    # import precalculated information about gaps in the alignment of a specific protein familiy
    gap_filename <- paste(path_to_folder_with_gap_files, tigr, ".js", sep="")
    gap <- read.table(gap_filename, header=F)
    
    # import precalculated dnds information of a specific protein family
    dnds_filename <- paste(path_to_folder_with_dnds_files, file, sep="")
    dnds <- read.table(dnds_filename, header=T)
    # sum up to get #N and #S
    dnds$dN <- dnds$leaf_nonsyn + dnds$tree_nonsyn
    dnds$dS <- dnds$leaf_syn + dnds$tree_syn
    dNall <-  dNall + sum(dnds$dN)
    dSall <-  dSall + sum(dnds$dS)
     
    # get sliding window information from whole sample 
    dN_win <- rollsum(dnds$dN, window_size,fill = list(NA, NULL, NA))
    dS_win <- rollsum(dnds$dS, window_size,fill = list(NA, NULL, NA))
    gap_win <- rollmean(gap$V1, window_size,fill = list(NA, NULL, NA))
    dnds_win <- dN_win/(dS_win+addsmallnumber)
    takethis_sample <- which(gap_win < gapthres)
    dnds_sample <- dnds_win[takethis_sample]
    dnds_sample_all <- c(dnds_sample_all,dnds_sample) 
<<<<<<< HEAD
    return_list <- list("N" = dNall, "S" = dSall, "dnds_window" = dnds_sample_all)
  } # for
  return(return_list)
} # function


=======
  } # for
  return(list(dNall, dSall, dnds_sample_all))
} # function


#### create sequence plot for a specific familiy ####
createSequencePlot <- function(tigr) {

  
  return(result)
}
>>>>>>> 79c713afd62f09ba7294774892c0778500191683
