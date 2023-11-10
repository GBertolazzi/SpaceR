
spatial_dependence_test <-  function(data,X, Y,  markers,alternative ="two.sided", B=5000, seed=NULL,  scatterplot=TRUE,  pch=16, colors=NULL,color_neg="blue", name_plot="scatterplot"){
  
  list.of.packages <- c( "TSdist",  "usedist") 
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){ install.packages(new.packages)}
  
  library(usedist)
  require(TSdist)

  lista = list() 
  n= numeric()
  
  for(i in 1:length(markers)){
    names(data)[which(  names(data)==markers[i])] = paste0("M", i)
    
    lista[[i]] = which(data[ , names(data)== paste0("M", i)]==1)
    
    n[i] =sum(data[ , names(data)== paste0("M", i)]) 
  }
  
  positive = paste0("M", c(1:length(markers)))
  names(lista) = positive
  #####
  
  ###  X coordinates
  if(length(X)==1){
    names(data)[which(  names(data)==X)] = "X" }
  
  if(length(X)==2){
    names(data)[which(  names(data)==X[1])] = "Xmin"
    names(data)[which(  names(data)==X[2])] = "Xmax"
    data$X = apply(cbind(data$Xmin, data$Xmax), 1, mean )   }
  
  ## Y coordinates
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
  
  ##
  coordinate = cbind(data$X, data$Y)
  
  ## 
  no_M1 = setdiff( c(1:nrow(data)), lista[[1]] )
  
  N = nrow(data) 
  
  n_M1 = length(lista[[1]])
  
  if(is.null(seed)==FALSE){ set.seed(seed)}
  
  d =dist(coordinate)
  
  A_all = as.matrix(d)
  
  ## 
  average_NN_observed = n_dp = numeric()
  
  for(j in 2:length(lista)){
    
    A = A_all[ lista[[1]] , lista[[j]] ]
    
    NN =numeric()
    
    for(i in 1:nrow(A)){
      NN[i] =    min(A[i,], na.rm = T)
    }
    
    ###
    average_NN_observed[j-1] =  mean(NN)
    
    n_dp[j] = length(which(NN==0)) 
    
    if(  n_dp[j]>0){
      cat(paste0("warning: there are ", n_dp[j]," ",  markers[1], "-", markers[j], " double positive cells \n"))}
    
  } # end cycle in j
  
  if( length(which(n_dp>0))>0 ){  
    cat("warning: the number of double positive cells will be fixed during the resampling procedure. \n \n")}
  
  ###############
  #### resampling
  B2 = matrix(NA, (length(lista)-1), 0)
  
  cat("Cell resampling... \n")
  for(k in 1:B){
    
    average_NN = numeric()
    for(j in 2:length(lista)){
      
      ## 
      R_double =  sample( lista[[1]], n_dp[j], replace = F )
      R_other =  sample( no_M1 , ( n[j] - n_dp[j] ), replace = F )
      
      A = A_all[ lista[[1]],   c(R_double, R_other )   ]
      
      NN = apply(A, 1, min)
      
      average_NN[j-1] =  mean(NN)
      
    } ## fine ciclo in j
    
    B2 = cbind(B2, average_NN)
    
    if(k == as.integer(B*20/100) ){  cat("20%  Done  \n")}
    if(k == as.integer(B*40/100) ){  cat("40%  Done  \n")}
    if(k == as.integer(B*60/100) ){  cat("60%  Done  \n") }
    if(k == as.integer(B*80/100) ){  cat("80%  Done  \n") }
    if(k == B ){  cat("100%  Done  \n") }  
    
  } ## end cycle in j
  
  cat("--- \n")
  ########
  ## 
  
  rand_average_NN =apply(B2, 1 , mean)
  
  R = round(average_NN_observed /rand_average_NN ,3)
  
  pval= numeric()
  
  if(alternative=="aggregation"){
    for(i in 1:(length(lista)-1)){
      pval[i] = length(which( B2[i, ]<=average_NN_observed[i] ))/ncol(B2) } }
  
  if(alternative =="segregation"){
    for(i in 1:(length(lista)-1)){
      pval[i] = length(which( B2[i, ]>=average_NN_observed[i] ))/ncol(B2) } }
  
  if(alternative =="two.sided"){
    B3 =   cbind(average_NN_observed , B2)
    Bz0 = apply(B3, 1, scale )
    
    Bz = t(Bz0[ -1, ])
    average_NN_observed_z = Bz0[1, ]
    
    for(i in 1:(length(lista)-1)){
      pval[i] = length(which( Bz[i, ]>=abs(average_NN_observed_z[i]) ))/ncol(B2) }
  } # 
  
  if(alternative !="two.sided" & alternative !="segregation" & alternative!="aggregation"){
    stop(paste0("unknown option: alternative = ",alternative ))   }
  
  pval_pr = ifelse(pval==0, 1/(B+1), pval)
  significativ =  ifelse(pval_pr<0.001, "***", 
                         ifelse(pval_pr<0.01, "**", 
                                ifelse(pval_pr<0.05, "*", 
                                       ifelse(pval_pr<0.1, ".", " "))))
  
  pval = ifelse(pval==0,paste0("<", (1/B)), pval )
  
  ris = data.frame( markers = markers[-1], n = n[-1],  average_NN= round( average_NN_observed,2),
                    rand_average_NN= round(rand_average_NN,2), R, pval  )
  
  
  ris_printing = data.frame(ris[, -1], significativ)
  rownames(ris_printing) = markers[-1]
  names(ris_printing)[ncol(ris_printing)] = " "
  
  cat( paste0("Average distances of nearest neighbor cells to ", markers[1]) , "(n = ", n[1],"):  \n")
  print(ris_printing )
  
  cat("--- \n")
  cat(paste0("alternative hypothesis: ",alternative ," \n" )) 
  cat("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 \n")
  
  if(scatterplot==TRUE){
    # scatterplot
    if(is.null(colors)==TRUE){ 
      colori=c("red", "green", "orange", "violet","yellow", "forestgreen") 
      colori = colori[1:length(lista)]
    }else{
      colori = colors }
    
    pdf(paste0(name_plot, ".pdf")) 
    
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
  
} #end function
