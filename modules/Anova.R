CalculateAnova <- function (transdf) {
  
  potato <- transdf %>% gather(key = m.z, value = intensity, -Group)
  
  anovaPTable <- data.frame()
  
  for (i in unique(potato$m.z)) {
    
    borg <- summary(aov(intensity ~ Group, data = potato[which(potato$m.z == i),]))
    glorp <- data.frame(m.z = i, AnovaPvalue = borg[[1]][["Pr(>F)"]][1])
    anovaPTable <- rbind(anovaPTable, glorp)
    
  }
  
  return(anovaPTable)
  
}