removeGroupMissingness <- function (workingDf, groupIdentities, wholeMissingness, groupMissingness) {

  workingDf <- workingDf[rowMeans(workingDf == 0) < (wholeMissingness/100), , drop = FALSE]

  tdf <- cbind(Group = groupIdentities$Group, as.data.frame(t(workingDf)))
  
  tdf <- tdf %>%     
    gather(key = m.z, value = intensity, -Group) %>%  #melt the data into a long format
    group_by(Group, m.z) %>% #group intensities by metabolite and cohort
    dplyr::summarise(perZ = sum(intensity==0)/length(intensity) < (groupMissingness/100)) %>%
    spread(Group, perZ)
  
  tdf <- tdf[rowSums(tdf[,2:length(tdf)]) > 0, 1]
  
  workingDf$m.z <- rownames(workingDf)
  workingDf <- merge(tdf, workingDf, by = "m.z", drop = FALSE)
  rownames(workingDf) <- workingDf$m.z
  workingDf$m.z <- NULL

  workingDf <- na.omit(workingDf)
  
  return(workingDf)
  
}


paretoScale <- function(df) {

  rowmean <- apply(df, 1, mean)
  rowsd <- apply(df, 1, sd)
  rowsqrtsd <- sqrt(rowsd)
  rv <- sweep(df, 1, rowmean, "-") #subtract row mean
  rv <- sweep(df, 1, rowsqrtsd, "/") #divide by square root of sd
  return(rv)
  
}

logTransform <- function(df) {

  df[df == 0] <- 1 #Impute zeros
  df <- log2(df) #Log transform
  return(df)
  
}

meanCenter <- function(df) {

  row_means <- rowMeans(df)
  
  # Subtract row means from each element in the corresponding row
  centered_df <- t(t(df) - row_means)
  
  # Convert the result back to a data frame
  centered_df <- as.data.frame(centered_df)
  
  return(centered_df)
}
