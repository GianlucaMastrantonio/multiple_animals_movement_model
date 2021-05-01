
### Parameter to set
DIR_CODE = ""

### Libraries
library(ggplot2)
library(xtable)
library(reshape2)
library(hrbrthemes)


#### #### #### #### #### #### #### ####
####  Libraries
#### #### #### #### #### #### #### ####

library(ggplot2)
#library(ggmosaic)
library(moveHMM)

### ### ### ### ### ###
### Functions
### ### ### ### ### ###

findmode = function(x){
        TT = table(as.vector(x))
        return(as.numeric(names(TT)[TT==max(TT)][1]))
}
dmnorm=function (x, mean = rep(0, d), varcov, log = FALSE)
{
    d <- if (is.matrix(varcov))
        ncol(varcov)
    else 1
    if (d > 1 & is.vector(x))
        x <- matrix(x, 1, d)
    n <- if (d == 1)
        length(x)
    else nrow(x)
    X <- t(matrix(x, nrow = n, ncol = d)) - mean
    Q <- apply((solve(varcov) %*% X) * X, 2, sum)
    logDet <- sum(logb(abs(diag(qr(varcov)[[1]]))))
    logPDF <- as.vector(Q + d * logb(2 * pi) + logDet)/(-2)
    if (log)
        logPDF
    else exp(logPDF)
}


#### #### #### #### #### #### ####
#### #### DIRECTORIES
#### #### #### #### #### #### ####
DIR_DATA = paste(DIR_CODE,"output/",sep="")
DIR_OUT  = paste(DIR_CODE,"output/Figures/",sep="")

### ### ### ### ### ### ### ### ###
### CODE
### ### ### ### ### ### ### ### ###
RDATA_NAME = "MODELOUT_Prop_3210_1_1_K=100.Rdata"

### color palette and ggplot objects
cbPalette     = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00")
cbPalette     = c(cbPalette,cbPalette)
GgplotTheme   = theme(
    axis.text.x   = element_text(face="bold",size=25),
    axis.text.y   = element_text(face="bold",size=25),
    axis.title.x  = element_text(face="bold",size=25),
    axis.title.y  = element_text(face="bold",size=25)
)

### ### LOAD ### ###
load(paste(DIR_DATA,RDATA_NAME,sep=""))


nanim           = length(MCMCMODEL_OUT)
vecOrdK         = list()
IndexMaxMCMC    = list()
Kmcmc           = list()
ZMAP            = list()

### FIgures - Posterior distribution of the number of behaviors - Not in the paper
for(ian in 1:nanim)
{
  ModelOUT = MCMCMODEL_OUT[[ian]]

  WWzeta = apply(ModelOUT$zeta,1,function(x) length(unique(x)))
  Lzeta  = length(unique(WWzeta))
  Uzeta  = unique(WWzeta)

  Data_ggplot = data.frame(K=WWzeta)

  p = ggplot(Data_ggplot, aes(x= K))+geom_bar(aes(y =..prop..),fill="steelblue")+ GgplotTheme+xlab("K")+ylab("Probability")+ylim(0,1)+scale_x_continuous( breaks = 1:10, limits=c(1,10))

  pdf(paste(DIR_OUT,"barplot_an",ian,".pdf",sep=""))
  print(p)
  dev.off()

  TT            = table(WWzeta)
  Kmcmc[[ian]]  = as.numeric(names(which.max(TT)))

  IndexMaxMCMC[[ian]] = which(WWzeta==Kmcmc[[ian]])

  ZMAP[[ian]] = apply(ModelOUT$zeta[IndexMaxMCMC[[ian]],],2,findmode)

  TT = table(ZMAP[[ian]])
  vecOrdK[[ian]] =  as.numeric(names(TT[order(-TT)]))
}


### Tables - Posterior estimates
for(ian in 1:nanim)
{
  ModelOUT = MCMCMODEL_OUT[[ian]]
  WM       =  vecOrdK[[ian]]

  print("")
  print("")
  print("")
  print(paste("Model",ian))
  for(i in 1:2)
  {
    dd = paste("$backslash mu_{j,",i,"}$ ",sep="")
    ff = "(CI) "
    for(k in WM)
    {
      dd = paste(dd, " & ", round(mean(ModelOUT$mu[,i,k]),3),sep = )
      ff = paste(ff, " & ","(",round(quantile(ModelOUT$mu[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$mu[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
    }
    dd = paste(dd, " \\")
    ff = paste(ff, " \\")
    print(dd)
    print(ff)
  }
  for(i in 1:2)
  {
    dd = paste("$backslash eta_{j,",i,"}$ ",sep="")
    ff = "(CI) "
    for(k in WM)
    {
      dd = paste(dd, " & ", round(mean(ModelOUT$eta[,i,k]),3),sep = )
      ff = paste(ff, " & ","(",round(quantile(ModelOUT$eta[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$eta[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
    }
    dd = paste(dd, " \\")
    ff = paste(ff, " \\")
    print(dd)
    print(ff)
  }
  for(i in 1:1)
  {
    dd = paste("$backslash nu_{j}$ ",sep="")
    ff = "(CI) "
    for(k in WM)
    {
      dd = paste(dd, " & ", round(mean(ModelOUT$nu[,i,k]),3),sep = )
      ff = paste(ff, " & ","(",round(quantile(ModelOUT$nu[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$nu[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
    }
    dd = paste(dd, " \\")
    ff = paste(ff, " \\")
    print(dd)
    print(ff)
  }

  for(i in 1:1)
  {
    dd = paste("$backslash rho_{j}$ ",sep="")
    ff = "(CI) "
    for(k in WM)
    {
      dd = paste(dd, " & ", round(mean(ModelOUT$rho[,i,k]),3),sep = )
      ff = paste(ff, " & ","(",round(quantile(ModelOUT$rho[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$rho[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
    }
    dd = paste(dd, " \\")
    ff = paste(ff, " \\")
    print(dd)
    print(ff)
  }


  for(i in 1:2)
  {
    for(j in i:2)
    {
      dd = paste("$backslash boldsymbol{backslash Sigma}_{", i,",",j, "}$ ",sep="")
      ff = "(CI) "
      for(k in WM)
      {
        dd = paste(dd, " & ", round(mean(ModelOUT$sigma[,(i-1)*2+j,k]),3),sep = )
        ff = paste(ff, " & ","(",round(quantile(ModelOUT$sigma[,(i-1)*2+j,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$sigma[,(i-1)*2+j,k], probs=c(1-0.025)),3)  ,   ")",sep="")
      }
      dd = paste(dd, " \\")
      ff = paste(ff, " \\")
      print(dd)
      print(ff)
    }


  }

  ik=1
  for(k in WM)
  {
    dd = paste("$backslash boldsymbol{backslash pi}_{", ik, "}$ ",sep="")
    ff = "(CI) "
    ii = 1
    for(i in WM)
    {
      dd = paste(dd, " & ", round(mean(ModelOUT$pi[,i,k]),3),sep = )
      ff = paste(ff, " & ","(",round(quantile(ModelOUT$pi[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$pi[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
      ii = ii+1
    }
    dd = paste(dd, " \\")
    ff = paste(ff, " \\")
    print(dd)
    print(ff)
    ik = ik+1
  }
  ik=1
  dd = paste("$n_{j}$ ",sep="")
  #ff = "(CI) "
  ii = 1
  for(i in WM)
  {
    dd = paste(dd, " & ", sum(ZMAP[[ian]]==i),sep = )
    #ff = paste(ff, " & ","(",round(quantile(ModelOUT$pi[,i,k], probs=c(0.025)),3), " ",round(quantile(ModelOUT$pi[,i,k], probs=c(1-0.025)),3)  ,   ")",sep="")
    ii = ii+1
  }
  dd = paste(dd, " \\")
  #ff = paste(ff, " \\")
  print(dd)
  #print(ff)
  ik = ik+1

}

### observed spatial location
for(ian in 1:nanim)
{
  ModelOUT = MCMCMODEL_OUT[[ian]]
  WM       =  vecOrdK[[ian]]

  ik = 1
  # for(i in WM)
  # {
  #   W = ZMAP[[ian]]==i
  #
  #   DataZ = data.frame(Longitude = DataCoords[W,(ian-1)*2+1], Latitude = DataCoords[W,(ian-1)*2+2] )
  #
  #   p = ggplot(DataZ, aes(Longitude, Latitude, shape=ik))
  #   p = p +geom_point(size = 2.5, shape=16,color=cbPalette[ik])
  #   p = p+theme(
  #     axis.text.x = element_text(face="bold",size=25),
  #     axis.text.y = element_text(face="bold",size=25),
  #     axis.title.x = element_text(face="bold",size=25),
  #     axis.title.y = element_text(face="bold",size=25),
  #     legend.text = element_blank(),
  #     legend.title = element_blank(),
  #     legend.position = "none"
  #   )+xlim(
  #     min(DataCoords[,seq(1,ncol(DataCoords),by=2)], na.rm=T),              max(DataCoords[,seq(1,ncol(DataCoords),by=2)], na.rm=T)
  #   ) +ylim(
  #     min(DataCoords[,seq(2,ncol(DataCoords),by=2)], na.rm=T), max(DataCoords[,seq(2,ncol(DataCoords),by=2)], na.rm=T)
  #   )
  #   p
  #   pdf(paste(DIR_OUT ,"DataPost_","Anim",ian,"k=",ik,".pdf",sep=""))
  #   print(p)
  #   dev.off()
  #   ik = ik+1
  # }
  DataZ = data.frame(Longitude = DataCoords[,(ian-1)*2+1], Latitude = DataCoords[,(ian-1)*2+2] )

  p = ggplot(DataZ, aes(Longitude, Latitude, shape=ik))
  p = p +geom_point(size = 2.5, shape=16)
  p = p+theme(
    axis.text.x = element_text(face="bold",size=25),
    axis.text.y = element_text(face="bold",size=25),
    axis.title.x = element_text(face="bold",size=25),
    axis.title.y = element_text(face="bold",size=25),
    legend.text = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )+xlim(
    min(DataCoords[,seq(1,ncol(DataCoords),by=2)], na.rm=T),              max(DataCoords[,seq(1,ncol(DataCoords),by=2)], na.rm=T)
  ) +ylim(
    min(DataCoords[,seq(2,ncol(DataCoords),by=2)], na.rm=T), max(DataCoords[,seq(2,ncol(DataCoords),by=2)], na.rm=T)
  )
  p
  pdf(paste(DIR_OUT ,"DataPost_","Anim",ian,".pdf",sep=""))
  print(p)
  dev.off()


}


### Attractive points
DataZ = data.frame(Longitude = c(DataCoords[,1],DataCoords[,3],DataCoords[,5],DataCoords[,7],DataCoords[,9],DataCoords[,11]), Latitude = c(DataCoords[,1+1],DataCoords[,3+1],DataCoords[,5+1],DataCoords[,7+1],DataCoords[,9+1],DataCoords[,11+1]) , Dog= rep(c("Woody","Sherlock","Alvin","Rosie","Bear","Lucy") , each= nrow(DataCoords)))

pointdata = data.frame(x= c(-0.072,0.575 ), y=c(0.263,-0.38  ), SpatialPoint= c("First","Second") )

p = ggplot(DataZ, aes(Longitude, Latitude))
  p = p +geom_point(aes(shape = Dog, color =Dog),
 size = 1.1)
  p = p+theme(
    axis.text.x = element_text(face="bold",size=25),
    axis.text.y = element_text(face="bold",size=25),
    axis.title.x = element_text(face="bold",size=25),
    axis.title.y = element_text(face="bold",size=25),legend.title = element_text( size = 24),
    legend.text = element_text( size = 20)
   # legend.text = element_blank(),
   # legend.title = element_blank(),
   # legend.position = "none"
  )+xlim(
    min(DataCoords[,seq(1,ncol(DataCoords),by=2)], na.rm=T),              max(DataCoords[,seq(1,ncol(DataCoords),by=2)], na.rm=T)
  ) +ylim(
    min(DataCoords[,seq(2,ncol(DataCoords),by=2)], na.rm=T), max(DataCoords[,seq(2,ncol(DataCoords),by=2)], na.rm=T)
  )+annotate("point",y = pointdata[,2], x=pointdata[,1], size=3.5)+
  annotate("text", y = pointdata[1,2]+0.3, x=pointdata[1,1], label = c(
  "First Attr. Point"
  ), size = 7
)+
 annotate("text", y = pointdata[2,2]-0.3, x=pointdata[2,1], label = c(
  "Second Attr. Point"
  ), size = 7
)


pdf(paste(DIR_OUT ,"DataPost_All",".pdf",sep=""))
print(p)
dev.off()


#
# for(ian in 1:nanim)
# {
#   ModelOUT = MCMCMODEL_OUT[[ian]]
#   WM       =  vecOrdK[[ian]]
#
#   ik = 1
#   for(i in WM)
#   {
#     W = ZMAP[[ian]]==i
#
#     DataZ = data.frame(Longitude = DataCoords[W,(ian-1)*2+1], Latitude = DataCoords[W,(ian-1)*2+2] )
#
#     p = ggplot(DataZ, aes(Longitude, Latitude, shape=ik))
#     p = p +geom_point(size = 2.5, shape=16,color=cbPalette[ik])
#     p = p+theme(
#       axis.text.x = element_text(face="bold",size=25),
#       axis.text.y = element_text(face="bold",size=25),
#       axis.title.x = element_text(face="bold",size=25),
#       axis.title.y = element_text(face="bold",size=25),
#       legend.text = element_blank(),
#       legend.title = element_blank(),
#       legend.position = "none"
#     )+xlim(
#       min(DataCoords[,seq(1,ncol(DataCoords),by=2)], na.rm=T),              max(DataCoords[,seq(1,ncol(DataCoords),by=2)], na.rm=T)
#     ) +ylim(
#       min(DataCoords[,seq(2,ncol(DataCoords),by=2)], na.rm=T), max(DataCoords[,seq(2,ncol(DataCoords),by=2)], na.rm=T)
#     )
#     p
#     pdf(paste(DIR_OUT ,"DataPost_","Anim",ian,"k=",ik,".pdf",sep=""))
#     print(p)
#     dev.off()
#     ik = ik+1
#   }
#   DataZ = data.frame(Longitude = DataCoords[,(ian-1)*2+1], Latitude = DataCoords[,(ian-1)*2+2] )
#
#   p = ggplot(DataZ, aes(Longitude, Latitude, shape=ik))
#   p = p +geom_point(size = 2.5, shape=16)
#   p = p+theme(
#     axis.text.x = element_text(face="bold",size=25),
#     axis.text.y = element_text(face="bold",size=25),
#     axis.title.x = element_text(face="bold",size=25),
#     axis.title.y = element_text(face="bold",size=25),
#     legend.text = element_blank(),
#     legend.title = element_blank(),
#     legend.position = "none"
#   )+xlim(
#     min(DataCoords[,seq(1,ncol(DataCoords),by=2)], na.rm=T),              max(DataCoords[,seq(1,ncol(DataCoords),by=2)], na.rm=T)
#   ) +ylim(
#     min(DataCoords[,seq(2,ncol(DataCoords),by=2)], na.rm=T), max(DataCoords[,seq(2,ncol(DataCoords),by=2)], na.rm=T)
#   )
#   p
#   pdf(paste(DIR_OUT ,"DataPost_","Anim",ian,".pdf",sep=""))
#   print(p)
#   dev.off()
#
#
# }


### SImil parameters
Kmcmc[[1]] = 3
Kmcmc[[2]] = 3
Kmcmc[[3]] = 4
Kmcmc[[4]] = 3
Kmcmc[[5]] = 3
Kmcmc[[6]] = 3

SimMatMu = matrix(NA, ncol=sum(unlist(Kmcmc)),nrow=sum(unlist(Kmcmc)))
SimMatEta = matrix(NA, ncol=sum(unlist(Kmcmc)),nrow=sum(unlist(Kmcmc)))
SimMatNu = matrix(NA, ncol=sum(unlist(Kmcmc)),nrow=sum(unlist(Kmcmc)))
SimMatRho = matrix(NA, ncol=sum(unlist(Kmcmc)),nrow=sum(unlist(Kmcmc)))
SimMatSigma = matrix(NA, ncol=sum(unlist(Kmcmc)),nrow=sum(unlist(Kmcmc)))



for(ian_1 in 1:nanim)
{

  ModelOUT_1 = MCMCMODEL_OUT[[ian_1]]
  WM_1         =  vecOrdK[[ian_1]]
  Wani_1 = c(1,cumsum(unlist(Kmcmc))+1)[ian_1]:cumsum(unlist(Kmcmc))[ian_1]

  for(ian_2 in 1:nanim)
  {
    ModelOUT_2 = MCMCMODEL_OUT[[ian_2]]
    WM_2        =  vecOrdK[[ian_2]]
    Wani_2 = c(1,cumsum(unlist(Kmcmc))+1)[ian_2]:cumsum(unlist(Kmcmc))[ian_2]


    #####
    Wsim = Reduce(intersect, list(IndexMaxMCMC[[ian_1]] ,IndexMaxMCMC[[ian_2]]))
    for(irow in 1:Kmcmc[[ian_1]])
    {

      k_lev1 = MCMCMODEL_OUT[[ian_1]]$z_mu[Wsim,vecOrdK[[ian_1]][irow]]
      Wrho0_1 = MCMCMODEL_OUT[[ian_1]]$rho[Wsim,1,vecOrdK[[ian_1]][irow]]==0
      Wrho1_1 = MCMCMODEL_OUT[[ian_1]]$rho[Wsim,1,vecOrdK[[ian_1]][irow]]==1
      k_lev1[Wrho1_1] = 2001
      for(icol in 1:Kmcmc[[ian_2]])
      {
        k_lev2 = MCMCMODEL_OUT[[ian_2]]$z_mu[Wsim,vecOrdK[[ian_2]][icol]]
        Wrho0_2 = MCMCMODEL_OUT[[ian_2]]$rho[Wsim,1,vecOrdK[[ian_2]][icol]]==0
        Wrho1_2 = MCMCMODEL_OUT[[ian_2]]$rho[Wsim,1,vecOrdK[[ian_2]][icol]]==1
        k_lev2[Wrho1_2] = 2002

        SimMatMu[Wani_1[irow], Wani_2[icol]] = mean(k_lev2==k_lev1)
      }
    }
    for(irow in 1:Kmcmc[[ian_1]])
    {
      k_lev1 = MCMCMODEL_OUT[[ian_1]]$z_eta[Wsim,vecOrdK[[ian_1]][irow]]
      Wrho0_1 = MCMCMODEL_OUT[[ian_1]]$rho[Wsim,1,vecOrdK[[ian_1]][irow]]==0
      Wrho1_1 = MCMCMODEL_OUT[[ian_1]]$rho[Wsim,1,vecOrdK[[ian_1]][irow]]==1
      k_lev1[Wrho0_1] = 2001
      for(icol in 1:Kmcmc[[ian_2]])
      {
        k_lev2 = MCMCMODEL_OUT[[ian_2]]$z_eta[Wsim,vecOrdK[[ian_2]][icol]]
        Wrho0_2 = MCMCMODEL_OUT[[ian_2]]$rho[Wsim,1,vecOrdK[[ian_2]][icol]]==0
        Wrho1_2 = MCMCMODEL_OUT[[ian_2]]$rho[Wsim,1,vecOrdK[[ian_2]][icol]]==1
        k_lev2[Wrho0_2] = 2002
        SimMatEta[Wani_1[irow], Wani_2[icol]] = mean(k_lev2==k_lev1)
      }
    }
    for(irow in 1:Kmcmc[[ian_1]])
    {
      k_lev1 = MCMCMODEL_OUT[[ian_1]]$z_nu[Wsim,vecOrdK[[ian_1]][irow]]
      Wrho0_1 = MCMCMODEL_OUT[[ian_1]]$rho[Wsim,1,vecOrdK[[ian_1]][irow]]==0
      Wrho1_1 = MCMCMODEL_OUT[[ian_1]]$rho[Wsim,1,vecOrdK[[ian_1]][irow]]==1
      k_lev1[Wrho1_1] = 2001
      for(icol in 1:Kmcmc[[ian_2]])
      {
        k_lev2 = MCMCMODEL_OUT[[ian_2]]$z_nu[Wsim,vecOrdK[[ian_2]][icol]]
        Wrho0_2 = MCMCMODEL_OUT[[ian_2]]$rho[Wsim,1,vecOrdK[[ian_2]][icol]]==0
        Wrho1_2 = MCMCMODEL_OUT[[ian_2]]$rho[Wsim,1,vecOrdK[[ian_2]][icol]]==1
        k_lev2[Wrho1_2] = 2002

        SimMatNu[Wani_1[irow], Wani_2[icol]] = mean(k_lev2==k_lev1)
      }
    }
    for(irow in 1:Kmcmc[[ian_1]])
    {
      k_lev1 = MCMCMODEL_OUT[[ian_1]]$z_rho[Wsim,vecOrdK[[ian_1]][irow]]
      for(icol in 1:Kmcmc[[ian_2]])
      {
        k_lev2 = MCMCMODEL_OUT[[ian_2]]$z_rho[Wsim,vecOrdK[[ian_2]][icol]]

        SimMatRho[Wani_1[irow], Wani_2[icol]] = mean(k_lev2==k_lev1)
      }
    }
    for(irow in 1:Kmcmc[[ian_1]])
    {
      k_lev1 = MCMCMODEL_OUT[[ian_1]]$z_sigma[Wsim,vecOrdK[[ian_1]][irow]]
      for(icol in 1:Kmcmc[[ian_2]])
      {
        k_lev2 = MCMCMODEL_OUT[[ian_2]]$z_sigma[Wsim,vecOrdK[[ian_2]][icol]]

        SimMatSigma[Wani_1[irow], Wani_2[icol]] = mean(k_lev2==k_lev1)
      }
    }


  }
}

NAMEPAR = c("mu","eta","rho","nu","sigma")

diag(SimMatEta) = NA
diag(SimMatSigma) = NA
diag(SimMatMu) = NA
diag(SimMatSigma) = NA
diag(SimMatNu) = NA
diag(SimMatRho) = NA
# Index = c(
# 	c(1,2,3),
# 	c(1,2,3)+length(vecOrdK[[1]]),
# 	c(1,2,3,4)+length(vecOrdK[[1]])+length(vecOrdK[[2]]),
# 	c(1,2,3)+length(vecOrdK[[1]])+length(vecOrdK[[2]])+length(vecOrdK[[3]]),
# 	c(1,2,3)+length(vecOrdK[[1]])+length(vecOrdK[[2]])+length(vecOrdK[[3]])+length(vecOrdK[[4]]),
# 	c(1,2,3)+20
#
# )

#image(SimMatMu)
# SimMatMu = SimMatMu[Index,Index]
# SimMatEta = SimMatEta[Index,Index]
# SimMatRho = SimMatRho[Index,Index]
# SimMatNu = SimMatNu[Index,Index]
# SimMatSigma = SimMatSigma[Index,Index]


# SimMatMu = SimMatMu[1:sum(unlist(Kmcmc)),1:sum(unlist(Kmcmc))]
# SimMatEta = SimMatEta[1:sum(unlist(Kmcmc)),1:sum(unlist(Kmcmc))]
# SimMatRho = SimMatRho[1:sum(unlist(Kmcmc)),1:sum(unlist(Kmcmc))]
# SimMatNu = SimMatNu[1:sum(unlist(Kmcmc)),1:sum(unlist(Kmcmc))]
# SimMatSigma = SimMatSigma[1:sum(unlist(Kmcmc)),1:sum(unlist(Kmcmc))]

for(i in 1:5)
{
  if(i==1){ddd = SimMatMu}
  if(i==2){ddd = SimMatEta}
  if(i==3){ddd = SimMatRho}
  if(i==4){ddd = SimMatNu}
  if(i==5){ddd = SimMatSigma}
  dataGG = melt(ddd)
  dataGG$Var2 = max(dataGG$Var2)-dataGG$Var2+1

  P = ggplot(dataGG , aes(Var1,Var2, fill=value)) + geom_raster()+scale_fill_gradient2(low="white",mid="orange", high="red", trans = "exp",limits=c(0,1), midpoint = exp(0.7))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    #axis.line = element_line(colour = "black"),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    legend.title = element_text( size = 24),
    legend.text = element_text( size = 20),
    axis.ticks.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
  )+labs(fill = "Prob.")+
  geom_hline(yintercept=1:(nrow(SimMatNu)+1)-0.5,size=0.1)+
  geom_vline(xintercept=(1:(nrow(SimMatNu)+1)-0.5) ,size=0.1)+
  geom_hline(yintercept=nrow(SimMatNu)+1-(c(1,cumsum(unlist(Kmcmc))+1)-0.5),size=1.1)+
  geom_vline(xintercept=((c(1,cumsum(unlist(Kmcmc))+1)-0.5)),size=1.1)+
  coord_fixed(ratio = 1)+annotate("text", x = 1:19, y=20.3, label = c(
  "B11","B12","B13",
  "B21","B22","B23", "B24",
  "B31","B32","B33",
  "B41","B42","B43",
  "B51","B52","B53",
  "B61","B62","B63"

  ),angle = -90, size = 5
)+annotate("text", y = rev(1:19), x=-0.3, label = c(
  "B11","B12","B13",
  "B21","B22","B23", "B24",
  "B31","B32","B33",
  "B41","B42","B43",
  "B51","B52","B53",
  "B61","B62","B63"

  ),angle = 0, size = 5
)+annotate("text", y = rev(c(2,5,8,11.5,15,18)), x=-2, label = c(
  "Woody", "Sherl.","Alvin","Rosie","Bear","Lucy"

  ),angle = 90, size = 5
)+annotate("text", x = (c(2,5,8.5,12,15,18)), y=20+2, label = c(
  "Woody", "Sherl.","Alvin","Rosie","Bear","Lucy"

  ),angle = 0, size = 5
)

  png(paste(DIR_OUT ,"Sim_",NAMEPAR[i],".png",sep=""))
  print(P)
  dev.off()
}




###
P_SL  = SimMatEta*SimMatSigma
P_OU  = SimMatMu*SimMatSigma*SimMatNu
P_TOT = P_SL
P_TOT[c(3,6,8,9,10,13,16,19),c(3,6,8,9,10,13,16,19)] = P_OU[c(3,6,8,9,10,13,16,19),c(3,6,8,9,10,13,16,19)]


LL = 1
dataGG = melt(P_TOT)
  dataGG$Var2 = max(dataGG$Var2)-dataGG$Var2+1

  P = ggplot(dataGG , aes(Var1,Var2, fill=value)) + geom_raster()+scale_fill_gradient2(low="white",mid="orange", high="red", trans = "exp",limits=c(0,  LL), midpoint = exp(LL*0.7))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    #axis.line = element_line(colour = "black"),
    axis.title.x=element_blank(),
    legend.title = element_text( size = 24),
    legend.text = element_text( size = 20),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
  )+labs(fill = "Prob.")+
  geom_hline(yintercept=1:(nrow(SimMatNu)+1)-0.5,size=0.1)+
  geom_vline(xintercept=(1:(nrow(SimMatNu)+1)-0.5) ,size=0.1)+
  geom_hline(yintercept=nrow(SimMatNu)+1-(c(1,cumsum(unlist(Kmcmc))+1)-0.5),size=1.1)+
  geom_vline(xintercept=((c(1,cumsum(unlist(Kmcmc))+1)-0.5)),size=1.1)+
  coord_fixed(ratio = 1)+annotate("text", x = 1:19, y=20.3, label = c(
  "B11","B12","B13",
  "B21","B22","B23", "B24",
  "B31","B32","B33",
  "B41","B42","B43",
  "B51","B52","B53",
  "B61","B62","B63"

  ),angle = -90, size = 5
)+annotate("text", y = rev(1:19), x=-0.3, label = c(
  "B11","B12","B13",
  "B21","B22","B23", "B24",
  "B31","B32","B33",
  "B41","B42","B43",
  "B51","B52","B53",
  "B61","B62","B63"

  ),angle = 0, size = 5
)+annotate("text", y = rev(c(2,5,8,11.5,15,18)), x=-2, label = c(
  "Woody", "Sherl.","Alvin","Rosie","Bear","Lucy"

  ),angle = 90, size = 5
)+annotate("text", x = (c(2,5,8.5,12,15,18)), y=20+2, label = c(
  "Woody", "Sherl.","Alvin","Rosie","Bear","Lucy"

  ),angle = 0, size = 5
)


  png(paste(DIR_OUT ,"Sim_Proc",".png",sep=""))
  print(P)
  dev.off()


ian = 1
nseq = 400
thetaseq = seq(-pi,pi,length.out=nseq+1)[-(nseq+1)]
rseq = seq(0,800,length.out=nseq)+0.00000000000000000000000000000001
diff1 = thetaseq[2]-thetaseq[1]
diff2 = rseq[2]-rseq[1]
Vecy  = matrix(NA, ncol=2, nrow=nseq*nseq)
Rvec     = rep(rseq, each=nseq)
Vecy[,1] = rep(rseq, each=nseq)*cos(rep(thetaseq,times=nseq))
Vecy[,2] = rep(rseq, each=nseq)*sin(rep(thetaseq,times=nseq))
ian = 1


for(ian in 1:nanim)
{

  ModelOUT = MCMCMODEL_OUT[[ian]]
  WM         =  vecOrdK[[ian]]

  sdsd = -200
  Densmat_t1 = list()
  Densmat_t2 = list()
  for(k in 1:length(vecOrdK[[ian]]))
  {
      Densmat_t1[[k]] = matrix(0, ncol=nseq,nseq)
      Densmat_t1[[k]] = matrix(0, ncol=nseq,nseq)
  }

  Wrho = apply(ModelOUT$rho[,1,]>0.9,2,mean)>0.9

  WWrho = which(Wrho[WM]==1)
  if(ian == 3)
  {
    WWrho = c(1)
  }
  if(ian == 4)
  {
    WWrho = c(1,2)
  }



  nmcmc          = dim(ModelOUT$mu)[1]
  if(length(WWrho)>0)
  {
    for(imcmc in 1:nmcmc)
    {
      k1 = 0
      for(k in WM[WWrho])
      {

        k1 = which(WM==k)
        mu = ModelOUT$eta[imcmc,,k]*SDtot*2*100000
        sigma = matrix(ModelOUT$sigma[imcmc,,k],2)*SDtot*2*100000*SDtot*2*100000
        Densmat_t1[[k1]][,] = Densmat_t1[[k1]][,]+matrix(dmnorm(Vecy,mu, sigma)*Rvec, ncol=nseq)/nmcmc

      }

    }
  }


  Densmat_2 = Densmat_t1



#lapply(lapply(Densmat_2,function(x){x*diff2} ),rowSums)
  DataTR = data.frame(
    Theta = c(unlist(lapply(lapply(Densmat_2,function(x){x*diff2} ),rowSums))),
    R     = c(unlist(lapply(lapply(Densmat_2,function(x){x*diff1} ),colSums))),
    Gruppo = as.factor(paste("B",ian,rep(1:length(vecOrdK[[ian]]),each =nseq ),sep="" ) ),
    thetaseq = rep(thetaseq, times=length(vecOrdK[[ian]])),
    rseq = rep(rseq, times=length(vecOrdK[[ian]]))
  )

  We = which(Wrho[vecOrdK[[ian]]]==1)
  if(ian == 4)
  {
    We = c(1,2)
  }
  DataTR  = DataTR[DataTR$Gruppo%in%as.character(paste("B",ian,We,sep="")),]
  #vecOrdK[[ian]]
#summary(DataTR)
  p = ggplot(DataTR,aes(x=thetaseq, y =Theta,group = Gruppo   ))+
  geom_line(aes(color=Gruppo),size=1.5)+ scale_fill_discrete(name="Bahavior")+theme(
    axis.text.y = element_text(face="bold",size=25),
    axis.text.x = element_text(face="bold",size=25),
    axis.title.x = element_text(face="bold",size=25),
    axis.title.y = element_text(face="bold",size=25),
    legend.text = element_text(face="bold",size=25),
    legend.title = element_blank()
  ) +ylim(0,0.23)+xlab("Turning-Angle")+ylab("Density") + labs(color="Behavior")+scale_color_manual(values=cbPalette)+guides(shape = guide_legend(override.aes = list(size = 5)))
   p
  #+scale_colour_brewer(
  #   type = "div",
  #   palette = "Spectral",
  #   direction = 1,
  #   aesthetics = "colour"
  #
  # )

  pdf(paste(DIR_OUT ,"an_",ian,"_Turn2.pdf",sep=""))
  print(p)
  dev.off()

  p = ggplot(DataTR,aes(x=rseq, y =R,group = Gruppo   ))+
  geom_line(aes(color=Gruppo),size=1.5)+ scale_fill_discrete(name="Bahavior")+theme(
    axis.text.y = element_text(face="bold",size=25),
    axis.text.x = element_text(face="bold",size=25),
    axis.title.x = element_text(face="bold",size=25),
    axis.title.y = element_text(face="bold",size=25),
    legend.text = element_text(face="bold",size=25),
    legend.title = element_blank()
  ) +ylim(0,max(DataTR$R))+xlab("Step-Length")+ylab("Density") + labs(color="Behavior")+scale_color_manual(values=cbPalette)+guides(shape = guide_legend(override.aes = list(size = 5)))
  p

  pdf(paste(DIR_OUT ,"an_",ian,"_Step2.pdf",sep=""))
  print(p)
  dev.off()

}
