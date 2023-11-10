
spatial_randomness_test <-  function(data,X, Y,  markers, B=5000, seed=NULL, scatterplot=TRUE,  pch=16, colors=NULL, color_neg="blue", name_plot="scatterplot"){
  
  list.of.packages <- c( "TSdist","usedist") 
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){ install.packages(new.packages)}
  
  library(usedist)
  library(TSdist)
  
  data$ID = paste0("c", c(1:nrow(data)))
  lista = list() 
  lunghezze = numeric()
  
  for(i in 1:length(markers)){
    names(data)[which(  names(data)==markers[i])] = paste0("M", i)
    
    lista[[i]] = which(data[ , names(data)== paste0("M", i)]==1)
    
    lunghezze[i] =sum(data[ , names(data)== paste0("M", i)]) 
  }
  
  positive = paste0("M", c(1:length(markers)))
  names(lista) = positive
  #####
  
  ### coordinate X
  if(length(X)==1){
    names(data)[which(  names(data)==X)] = "X" }
  
  if(length(X)==2){
    names(data)[which(  names(data)==X[1])] = "Xmin"
    names(data)[which(  names(data)==X[2])] = "Xmax"
    data$X = apply(cbind(data$Xmin, data$Xmax), 1, mean )   }
  
  ##coordinate Y
  if(length(Y)==1){
    names(data)[which(  names(data)==Y)] = "Y" }
  
  if(length(Y)==2){
    names(data)[which(  names(data)==Y[1])] = "Ymin"
    names(data)[which(  names(data)==Y[2])] = "Ymax"
    data$Y = apply(cbind(data$Ymin, data$Ymax), 1, mean )   }
  
  max_X = max(data$X)
  max_Y = max(data$Y)
  min_X = min(data$X)
  min_Y = min(data$Y)
  
  ##densita
  nrow(data)/((max_X-min_X)*(max_Y-min_Y))
  
  coordinate = cbind(data$X, data$Y)
  rownames(coordinate) = data$ID
  
  
  ## esclusione di marcatori che non presentano positività
  esclusion = which(lunghezze<2 )
  
  if(length(esclusion)>0){
    if(length(esclusion)==1){
      print(paste0(length(esclusion), " marker has been excluded"))
    }else{ print(paste0(length(esclusion), " markers have been excluded") )}
    lista = lista[-esclusion]
    positive = positive[-esclusion]
    markers = markers[-esclusion]
  }
  ###
  
  average_NN_observe = n= numeric()
  
  if(is.null(seed)==FALSE){ set.seed(seed)}
  
  for(j in 1:length(lista)){
    
    origin_row <- lista[[j]]
    
    sub_coord = coordinate[   origin_row ,  ]
    
    d =dist(sub_coord)
    
    n[j]  = length(origin_row)
    
    A = as.matrix(d)
    
    NN =numeric()
    
    for(i in 1:nrow(A)){
      A[i , which(colnames(A)==rownames(A)[i]  ) ] = NA
      NN[i] =    min(A[i,], na.rm = T)
    }
    
    ### average distance of the first neighbor in group j
    average_NN_observe[j] =  mean(NN)
    
  }
  
  ###############
  ####
  B2 = matrix(NA, length(lista), 0)
  
  cat("Spatial randomness testing... \n")
  
  for(k in 1:B){
    
    ### resampling
    lista_random = list()
    
    for(r in 1:length(lista)){
      lista_random[[r]] = sample(c(1:nrow(data)), length(lista[[r]]), replace = F )
    }
    
    names(lista_random)=  paste0("R_", names(lista))
    
    ## average distances in random samples
    average_NN = numeric()
    
    for(j in 1:length(lista_random)){
      
      origin_row <- lista_random[[j]]
      
      sub_coord = coordinate[  origin_row ,  ]
      
      d =dist(sub_coord)
      
      A = as.matrix(d)
      
      NN =numeric()
      
      for(i in 1:nrow(A)){
        A[i , which(colnames(A)==rownames(A)[i]  ) ] = NA
        NN[i] =    min(A[i,], na.rm = T)
      }
      
      average_NN[j] =  mean(NN)
      
    } ## end cycle in j
    
    B2 = cbind(B2, average_NN)
    
    if(k == as.integer(B*20/100) ){  cat("20%  Done  \n")}
    if(k == as.integer(B*40/100) ){  cat("40%  Done  \n")}
    if(k == as.integer(B*60/100) ){  cat("60%  Done  \n") }
    if(k == as.integer(B*80/100) ){  cat("80%  Done  \n") }
    if(k == B ){  cat("100%  Done  \n \n") }    
  } ## end cycle in k
  
  ## empirical measures
  
  rand_average_NN =apply(B2, 1 , mean)
  
  R = round(average_NN_observe /rand_average_NN ,3)
  
  pval= numeric()
  for(i in 1:length(lista)){
    pval[i] = length(which( B2[i, ]<=average_NN_observe[i] ))/ncol(B2) }
  
  pval_pr = ifelse(pval==0, 1/(B+1), pval)
  significativ =  ifelse(pval_pr<0.001, "***", 
                         ifelse(pval_pr<0.01, "**", 
                                ifelse(pval_pr<0.05, "*", 
                                       ifelse(pval_pr<0.1, ".", " "))))
  
  pval = ifelse(pval==0,paste0("<", (1/B)), pval )
  
  ris = data.frame(marker= markers,n,  average_NN= round( average_NN_observe,2),
                    rand_average_NN= round(rand_average_NN,2), R, pval  )
  
  ris_printing = data.frame(ris[, -1], significativ)
  rownames(ris_printing) = markers
  names(ris_printing)[6] = " "
  
  print(ris_printing )
  
  cat("--- \n")
  cat("alternative hypothesis: spatial clustering")
  cat("\n")
  cat("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  
  
  if(scatterplot==TRUE){
    # scatterplot
    if(is.null(colors)==TRUE){ 
      colori=c("red", "green", "orange", "violet","yellow", "forestgreen") 
      colori = colori[1:length(lista)]
    }else{
      colori = colors }
    
    pdf(paste0(name_plot,".pdf")) 
    
    plot(data$X[1],data$Y[1],col="white",  xlim =c(min_X, max_X), ylim =c(min_Y, max_Y), 
         xlab = "X", ylab="Y",main="cell spatial distribution", pch=16 )
    
    for(j in 1:length(lista)){
      dataM = data[ lista[[j]] ,c("X", "Y")]
      points(dataM$X,dataM$Y, col=colori[j], pch=16)
    }
    
    negative = setdiff(c(1:nrow(data)), unlist(lista) )
    
    if(length(pch)==1){ 
      pch2= rep(pch, length(lista))
    }else{ pch=pch2 }
    
    if(length(negative)>0){
      dataM = data[negative, c("X", "Y")]
      points(dataM$X,dataM$Y, pch=1, col=color_neg) 
      legend("topleft", bg="white", legend= c(markers,"negative") ,
             col=c(colori, color_neg),pch = c(pch2, 1),  cex=1)
    }else{
      legend("topleft", bg="white", legend= markers ,
             col=colori ,pch = pch2,  cex=1)
    }
    
    dev.off()
  } # end plot
  
  
  return(ris)
  
} 
##################

