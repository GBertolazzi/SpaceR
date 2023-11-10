
distribution_equality_test = function(data, X , Y,   M1 , M2 , colM1="green" ,  colM2="red" , rm.double=FALSE, B=5000 , seed=NULL ){
  
  list.of.packages <- c( "TSdist",  "spatstat", "graphics") 
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){ install.packages(new.packages)}
  
#  require(ape)
  require(TSdist)
  require(spatstat)
  require(graphics)
  
  names(data)[which(  names(data)==M1)] = "M1"
  names(data)[which(  names(data)==M2)] = "M2"
  
  ## stop se ci sono doppie positive e rm.double=F
n_dp =   length(which(data$M1==1 & data$M2==1))
  
if(n_dp>0 & rm.double==FALSE ){
    stop("stop: double positive cells are not allowed")   }

  if(n_dp>0 & rm.double==TRUE ){
    cat("warning: ",n_dp, " double positive cells have been removed \n \n " )   }
  
  ## elimino le cellule double positive e doppie negative
  data = data[ -which((data$M1==1 & data$M2==1) | (data$M1==0 & data$M2==0)   ) , ]

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
  
  ## creazione della finestra
  owinf_obj = owin(xrange=c(min_X,max_X), yrange=c(min_Y,max_Y))
  
  data_M1 = subset(data, M2==0 & M1==1 )
  data_M2 = subset(data, M2==1 & M1==0 )
  
  n_M1 = nrow(data_M1) 
  n_M2 = nrow(data_M2)
  
  ## prob di M1
  prob = nrow(data_M1)/nrow(data) 
  
  pp_all = ppp(data$X, data$Y   , window = owinf_obj ) 
  pp_M1 = ppp(data_M1$X, data_M1$Y   , window = owinf_obj ) 
  pp_M2 = ppp(data_M2$X, data_M2$Y   , window = owinf_obj ) 
   
  count_all =  as.matrix(quadratcount(pp_all))
  count_M1 =  as.matrix(quadratcount(pp_M1))
  count_M2 =  as.matrix(quadratcount(pp_M2))
  
  # CD3 attesi
  E = (count_all*prob)- count_all*(1-prob)
  
  O = count_M1-count_M2
  
  chi =   sum(((O-E)^2)/abs(E), na.rm=T)
  
  R = matrix( paste0("obs = " ,O, " \n exp = ", round(E, 2)), nrow(E), ncol(E), byrow = F)
  R =  ifelse(R=="obs = 0 \n exp = 0", "", R)
  
  #### ciclo in k per ricavare la distribuzione nulla!
  chi_rand = numeric()
  
  data_rand = data

  if(is.null(seed)==FALSE){ set.seed(seed)}
  
  cat("Distribution equality testing... \n")
  
  for(k in 1:B){
    
    data_R = data[, c("M1", "M2")]
    data_R = data_R[ sample(c(1:nrow(data_R)), replace=F)  , ]
    
    data_rand$M1 = data_R$M1
    data_rand$M2 = data_R$M2
    
    data_M1_rand = data_rand[ which(data_rand$M2==0 & data_rand$M1==1 ) , ]
    data_M2_rand = data_rand[ which(data_rand$M2==1 & data_rand$M1==0 ) , ]
    
    pp_M1 = ppp(data_M1_rand$X, data_M1_rand$Y   , window = owinf_obj ) 
    pp_M2 = ppp(data_M2_rand$X, data_M2_rand$Y   , window = owinf_obj ) 
    
    count_M1 =  as.matrix(quadratcount(pp_M1))
    count_M2 =  as.matrix(quadratcount(pp_M2))
    
    O = count_M1-count_M2
    
    chi_rand[k] =   sum(((O-E)^2)/abs(E), na.rm=T)
    
    if(k == as.integer(B*20/100) ){  cat("20%  Done  \n")}
    if(k == as.integer(B*40/100) ){  cat("40%  Done  \n")}
    if(k == as.integer(B*60/100) ){  cat("60%  Done  \n") }
    if(k == as.integer(B*80/100) ){  cat("80%  Done  \n") }
    if(k == B ){  cat("100%  Done  \n \n") }
  } # fine del ciclo in k
  
  pvalue = length(which(chi_rand>=chi))/length(chi_rand)
  pvalue_print =   ifelse(pvalue==0, paste0("< ",1/B ) , paste0("= ",signif(pvalue,2)))
  
  index =  chi/mean(chi_rand) #new index
  
  comparison = paste(M1, M2, sep="-") 
  
  ris_lista = list(   n_M1=n_M1, n_M2=n_M2, M1_proportion = prob,Chi_statistic =chi, pvalue =pvalue, index=index)
  
    pdf(paste0("scatterplot_",comparison, ".pdf"), width = 15) 
  
  par(mfrow=c(1, 2))
  
  # scatterplot
  plot(data_M1$X,data_M1$Y, col=colM1, xlim =c(min_X, max_X), ylim =c(min_Y, max_Y), 
       xlab = "X", ylab="Y",main="cell spatial distribution", pch=16 )
  points(data_M2$X,data_M2$Y, col=colM2, pch=17)
  
  Rxall = seq(min_X,max_X, length.out = (ncol(R)+1))
  Ryall  = seq(min_Y,max_Y, length.out = (nrow(R)+1) )
  
  abline(v=Rxall, lty=3, lwd=4, col="grey")
  abline(h= Ryall, lty=3, lwd=4, col="grey")
  
  legend("topleft", bg="white", legend=c(M1, M2 , paste0("p-value ", pvalue_print) ) ,
         col=c(colM1, colM2, "white"),pch = c(16,17, 1),  cex=1)
  
  ## griglia attesi/osservati
  plot(data_M1$X[1:2],data_M1$Y[1:2], col="white", xlim =c(min_X, max_X),
       ylim =c(min_Y, max_Y), xlab = "X", ylab="Y", pch=16, main = "Observed and expected differences", axes = TRUE )
  
  Rxall = seq(min_X,max_X, length.out = (ncol(R)+1) )
  Ryall  = seq(min_Y,max_Y, length.out = (nrow(R)+1) )
  
  abline(v=Rxall, lty=3, lwd=4, col="grey")
  abline(h= Ryall, lty=3, lwd=4, col="grey")
  
  Rx = seq(min_X,max_X, length.out = (ncol(R)*2+1) )[c(2,4,6,8,10)]
  Ry   = seq(max_Y,min_Y, length.out = (nrow(R)*2+1) )[c(2,4,6,8,10)]
  
  for(i in 1:length(Ry)){
    for(j in 1:length(Rx)){
      text(y= Ry[i], x=Rx[j] , R[i,j] , cex=1)
    }
  }
  dev.off()

  cat("test.statistic = ",round(chi,2), ", p-value ", pvalue_print , " \n")
  cat("alternative hypothesis: ",M1, " and ", M2, "cells have different spatial distributions \n")
  cat("dissimilarity index = ",round(index,2) , " \n")
  
  ## risultati
   return(ris_lista)
  
} # fine funzione
