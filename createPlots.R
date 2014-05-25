
createSequencePlot <- function(tigr,global_N,global_S,global_dnds_window) {
  ### create sequence plot for a specific familiy ###  
  # requires package zoo, ggplot2, intervals 
  # 

  cat(paste("Drawing",tigr,"\n"))
  # path to precalculated files
  dnds_filename <- paste(dnds_folder,tigr,dnds_file_ending,sep="")
  js_filename <- paste(js_folder,tigr,js_file_ending,sep="")
  gap_filename <- paste(gap_folder,tigr,gap_file_ending,sep="")
  membrane_filename <- paste(membrane_folder,tigr,membrane_file_ending,sep="")
 
  # read files
  gap <- read.table(gap_filename,header=F)
  dnds <- read.table(dnds_filename,header=T)
  js <- read.table(js_filename,header=F)
  membrane <- read.table(membrane_filename,header=F,sep=";")
  try(rep <- read.table(paste(repeat_folder,tigr,repeat_file_ending,sep="")),silent=TRUE)
  
  if (!is.null(rep)){
    repmat <- rep(FALSE, nrow(gap)) 
    for (x in 1:nrow(rep)){
      try(repmat[rep$V1[x]:rep$V2[x]] <- "TRUE",TRUE)
    }
  }  # if there are repeats, create matrix for plotting
  if(is.null(repmat)){
    repmat$V1[1] <- 0
    repmat$V2[1] <- 0
    rep <- as.data.frame(repmat)
  }  # if there are no repeats, creating empty matrix to prevent plotting problems
  if(is.null(membrane)){
    membrane$V1[1] <- 0
    membrane$V2[1] <- 0
    membrane$V3[1] <- 0
  } # if there is no membrane prediction available, create empty matrix 
  
  # calculate #N/#S
  dnds$dN <- dnds$leaf_nonsyn + dnds$tree_nonsyn
  dnds$dS <- dnds$leaf_syn + dnds$tree_syn
  dnds$dnds <- dnds$dN / (dnds$dS + addsmallnumber)
  dnds$js <- js$V1
  dnds$gap <- gap$V1  
  
  # count how many rep are in sequence in %
  rep_percent <- length(which(repmat==TRUE))/length(repmat)
  
  #smoothing with rollmean
  dnds_smooth <-rollmean(dnds$dnds,3,fill="extend") # 2nd parameter defines smoothness 
  dn_smooth <-rollmean(dnds$dN,3,fill="extend") # 2nd parameter defines smoothness
  ds_smooth <-rollmean(dnds$dS,3,fill="extend") # 2nd parameter defines smoothness 
  js_smooth <-rollmean(dnds$js,3,fill="extend") # 2nd parameter defines smoothness 
  
  #create data frame without deletion of positions with many gaps 
  df <- data.frame(dnds$pos,dnds_smooth)
  df_true <- data.frame(dnds$pos,dnds$dnds)
  colnames(df) <- c( 'Position', 'dnds')
  colnames(df_true) <- c( 'Position', 'dnds')
  df_true.long <- melt(df_true,id="Position")
  df.long <- melt(df, id="Position") 
  
  df_true.long$roll_dnds <- NULL
  roll_dn <- rollsum(dnds$dN, window_size, fill = list(NA, NULL, NA),align = "left")
  roll_ds <- rollsum(dnds$dS, window_size, fill = list(NA, NULL, NA),align = "left")
  roll_dnds <- roll_dn / (roll_ds+0.5) 
  roll_gap <- rollmean(dnds$gap, window_size,fill = list(NA, NULL, NA),align = "left")
  df_true.long$roll_dnds <- coredata(roll_dnds)  
  df_true.long$roll_gap <- coredata(roll_gap) 
  df_true.long$roll_dn <- coredata(roll_dn)
  df_true.long$roll_ds <- coredata(roll_ds) 
  
  # perform Fisher test in sliding window 
  df_true.long$window_sig <- 0
  iter <- 1
  for (dnds_value in df_true.long$roll_dnds){
    if (!is.na(dnds_value)){
      if (!is.na(df_true.long$roll_gap[iter])){
        if (df_true.long$roll_gap[iter]< gapthres){
          ctable <-  matrix(c(df_true.long$roll_dn[iter],global_N,df_true.long$roll_ds[iter],global_S),nrow=2,dimnames=list(c("all","sample"),c("dn","ds"))) # 2*2 table
          pval <- fisher.test(ctable,alternative=c("greater"))[["p.value"]]
          pval.corrected <- p.adjust(pval, method="fdr",n = length(global_dnds_window))
          x <- ifelse(pval < 0.01,1,0) 
          if (x == 1){ 
            df_true.long$window_sig[iter:(iter+window_size-1)] <- df_true.long$window_sig[iter:(iter+window_size-1)] + 1  
          } # if this is significant
        }
      }
    } 
    iter <- iter + 1
  } # iterate over sliding window and perform a fisher test for each window with fdr testing correction
  
  df_true.long$window_sig[which(df_true.long$window_sig > 0)] <- 1
  
  # create df for the conservation plot
  df_cons <-  data.frame(dnds$pos,js_smooth)
  colnames(df_cons) <- c('Position','value')
  long2gap <-  match(df_true.long$Position,row.names(gap))
  df_true.long$gap <- gap[long2gap,]
  meandnds <- sum(dnds$dN) / sum(dnds$dS)
  
  # crate tables
  dnds_gap_corrected <- mean(dnds$dnds)
  mean_js <- mean(dnds$js)
  mean_rep <- rep_percent
  mean_gap <- mean(dnds$gap)
  
  # get the start end end pos from 1
  it2 <- 1
  k <- 1
  area_sig<- NULL
  ind <- 1
  for (k in 1:length(df_true.long$window_sig)){  
    if (df_true.long$window_sig[k] == 1) {  # this elem is 1
      if (k ==1){area_sig$start[ind] <- k} else{
        if (df_true.long$window_sig[k-1] == 0){ # the last elem not
          area_sig$start[ind] <- k
        }}
    }
    if (!k ==1){ # first element can't be the end
      if (df_true.long$window_sig[k] == 0) {  # this elem is 0 
        if (df_true.long$window_sig[k-1] == 1){ # the last elem not
          area_sig$end[ind] <- k-1
          ind <- ind + 1
        }
      } 
    }
    # check if last element is 1 and make an end 
    if (k == length(df_true.long$window_sig)){
      if(df_true.long$window_sig[k] == 1){
        area_sig$end[ind] <- k
      }
    }
    k <- k +1 
  }
  # to prevent a if else statement on the plotting part if there is no significant region
  if(is.null(area_sig)){
    area_sig$start[1] <- 0
    area_sig$end[1] <- 0
  } # if no significant cluster is found
  area_sig <- as.data.frame(area_sig)
  area_sig$col <- "significant"
  
  # create talbe with summary information of clusters
  expo <- data.frame(matrix(vector(), length(dnds_list), 12, dimnames=list(c(), c("Name", "cluster_num", "cluster_size","inside","outside","helix","inside_size","outside_size","helix_size","inside_sig_size","outside_sig_size","helix_sig_size"))), stringsAsFactors=F)
  expo$inside <- 0; expo$outside <- 0; expo$helix<- 0
  expo$inside_size <- 0; expo$outside_size <- 0; expo$helix_size <- 0
  expo$inside_sig_size <- 0; expo$outside_sig_size <- 0; expo$helix_sig_size <- 0; tiNum <- 1
  expo$Name[tiNum] <- tigr
  expo$cluster_num[tiNum] <- nrow(area_sig)
  expo$cluster_size[tiNum] <- length(which(df_true.long$window_sig==1)) / (nrow(df_true.long))
  if (expo$cluster_size[tiNum] == 0){expo$cluster_num[tiNum] <- 0} # correct the issue if there is a significant cluster on 0-0
  
  # seperate membrane prediction into the three categories 
  membrane_outside <- membrane[which(membrane$V1 == "outside"),]
  membrane_inside <- membrane[which(membrane$V1 == "inside"),]
  membrane_helix <- membrane[which(membrane$V1 == "TMhelix"),]  
  if(nrow(membrane_outside)> 0){expo$outside_size[tiNum] <- sum(membrane_outside$V3-membrane_outside$V2)}
  if(nrow(membrane_inside)> 0){expo$inside_size[tiNum] <- sum(membrane_inside$V3-membrane_inside$V2)}
  if(nrow(membrane_helix)> 0){expo$helix_size[tiNum] <- sum(membrane_helix$V3-membrane_helix$V2)}

  # check if clusters are in TMhelix or inside/outside
  # iterate over significant position
  windowInt <- Intervals(area_sig[,1:2])
  membraneInt <- Intervals(membrane[,2:3])
  overlap <- interval_overlap(windowInt,membraneInt)
  overlap_vec <- c(do.call("cbind", overlap)) # convert list into vector 
  coolOnes <- membrane[overlap_vec,]
  
  if (nrow(coolOnes) > 0){
    # iterate through overlapping elements and count++
    for (overlap_element in 1 : nrow(coolOnes)){
      if (coolOnes$V1[overlap_element] == "inside"){expo$inside[tiNum] <- expo$inside[tiNum] + 1} #inside++
      if (coolOnes$V1[overlap_element] == "outside"){expo$outside[tiNum] <- expo$outside[tiNum] + 1} #inside++
      if (coolOnes$V1[overlap_element] == "TMhelix"){expo$helix[tiNum] <- expo$helix[tiNum] + 1} #inside++  
    }
  }
  
  
  
 # return()
}