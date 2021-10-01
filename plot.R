
### Parameter to set
DIR_CODE = ""

# name of the .Rdata with the model's results
RDATA_NAME = ".Rdata"
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


library(ggplot2)
library(ggmosaic)
library(moveHMM)
library(car)
library(mvtnorm)
### Libraries
library(ggplot2)
library(xtable)
library(reshape2)
library(hrbrthemes)

library(ggplot2)
library(moveHMM)
library(mcclust.ext)
### ### ### ### ### ###
### Colors palette
### ### ### ### ### ###
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette <- rev(c(
  "#ffffbf",
"#d7191c",
"#fdae61",
"#abdda4",
"#2b83ba")
)

findmode = function(x){
                    TT = table(as.vector(x))
                    return(as.numeric(names(TT)[TT==max(TT)][1]))
            }

#### #### #### #### #### #### #### ####
#### Coordinate Systems
#### #### #### #### #### #### #### ####

angle = -1/8*pi
RotMat = matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)),nrow=2)
#RotMat = diag(1,2)
CoordAx1 =  c(-1,0)%*%RotMat
CoordAx2 =  c(1,0)%*%RotMat
CoordAy1 =  c(0,1.)%*%RotMat
CoordAy2 =  c(0,-0.5)%*%RotMat

PT = qplot(c(-10),10,geom="line", xlab="",ylab="", ylim=c(-1,1), xlim=c(-1,1))+theme(axis.text=element_text(size=25),
        axis.title=element_text(size=28,face="bold"))+scale_colour_continuous(guide = FALSE)+theme(legend.position="none")+ coord_fixed(ratio=1)

PT1 = PT + geom_segment(aes(x = CoordAx1[1], y = CoordAx1[2], xend = CoordAx2[1], yend = CoordAx2[2]), linetype =  2)+geom_segment(aes(x = CoordAy1[1], y = CoordAy1[2], xend = CoordAy2[1], yend = CoordAy2[2]), linetype =  2)
# +
# geom_segment(aes(x = 0, y = 0, xend = 0, yend = 1), linetype =  2)+
# geom_segment(aes(x = 0, y = 0, xend = 1, yend = 0), linetype =  2)

OBS = matrix(c(-0.8,0,0.8,0,0,0.6),nrow=3,ncol=2)%*%RotMat

PT2 = PT1+geom_line(aes(x=OBS[,1],y=OBS[,2]))+geom_point(aes(x=OBS[,1],y=OBS[,2],size=1))
PT2
Desx1 = rbind(c(0.8,0.6),cbind(0.8,0))%*%RotMat
Desx2 = rbind(c(0.8,0.6),cbind(0,0.6))%*%RotMat
Desx1_v2 = Desx1
Desx1_v2[2,2] = 0
Desx1_v2[2,1] = Desx1_v2[1,1]
Desx2_v2 = Desx2
Desx2_v2[2,1] = 0
Desx2_v2[2,2] = Desx2_v2[1,2]

PT3 = PT2+geom_line(aes(x=Desx1[,1],y=Desx1[,2]), linetype = 3)+geom_line(aes(x=Desx2[,1],y=Desx2[,2]), linetype = 3)
#+geom_line(aes(x=Desx1_v2[,1],y=Desx1_v2[,2]), linetype = 3) + geom_line(aes(x=Desx2_v2[,1],y=Desx2_v2[,2]), linetype = 3)


Xc 		= rbind(c(0.8,0.02),c(0.8,-0.02))
Yc 		= rbind(c(-0.02,0.6),c(0.02,0.6))
Xc 		= Xc%*%RotMat
Yc 		= Yc%*%RotMat
Xc_v2 		= rbind(c(0.8,0.02),c(0.8,-0.02))
Xc_v2[1,]    =  OBS[3,]
Xc_v2[1,1]    =  -0.02
Xc_v2[2,1]    =  0.02
Xc_v2[2,2]  = Xc_v2[2,2]=  OBS[3,2]

Yc_v2 		= rbind(c(0.8,0.02),c(0.8,-0.02))
Yc_v2[1,]    =  OBS[3,]
Yc_v2[1,2]    =  -0.02
Yc_v2[2,2]    =  0.02
Yc_v2[2,1]  = Yc_v2[2,1]=  OBS[3,1]

PT4 = PT3+geom_line(aes(x=Xc[,1],y=Xc[,2]), linetype = 1,size=1)+geom_line(aes(x=Yc[,1],y=Yc[,2]), linetype = 1,size=1)
#+geom_line(aes(x=Xc_v2[,1],y=Xc_v2[,2]), linetype = 1,size=1)+geom_line(aes(x=Yc_v2[,1],y=Yc_v2[,2]), linetype = 1,size=1)


PT4

PT5 =PT4+
geom_text(aes(x=Xc[2,1]+0.1,y=Xc[2,2]-0.06), label=deparse(bquote(paste('y'['t'['i']*','][1]))),parse=TRUE,size=7)+
geom_text(aes(x=Yc[2,1]-0.11,y=Yc[2,2]-0.04), label=deparse(bquote(paste('y'['t'['i']*','][2]))),parse=TRUE,size=7)
PT5




PT6 =PT5+
geom_text(aes(x=OBS[1,1]-0.1,y=OBS[1,2]+0.04), label=deparse(bquote(paste('s'['t'['i-1']]))),parse=TRUE,size=7)+
geom_text(aes(x=OBS[2,1]-0.12,y=OBS[2,2]+0.02), label=deparse(bquote(paste('s'['t'['i']]))),parse=TRUE,size=7)+
geom_text(aes(x=OBS[3,1]+0.12,y=OBS[3,2]-0.02), label=deparse(bquote(paste('s'['t'['i+1']]))),parse=TRUE,size=7)
PT6


cangle = seq(0.38,pi/2-0.5,length.out=100)
S1 = sin(cangle)*0.2
C1 = cos(cangle)*0.2

PT7 = PT6 +geom_line(aes(x=C1,y=S1), linetype = 1,size=0.4)+geom_text(aes(x=0.2,y=0.18), label=deparse(bquote(paste(theta['t'['i']]))),parse=TRUE,size=7)+
geom_text(aes(x=0.37,y=0.45), label=deparse(bquote(paste('r'['t'['i']]))),parse=TRUE,size=7)
PT7


PT8 = PT7+geom_segment(aes(x = OBS[1,1], y = OBS[1,2], xend = 1, yend = OBS[1,2]), linetype =  3)+
geom_segment(aes(x = OBS[2,1], y = OBS[2,2], xend = 1, yend = OBS[2,2]), linetype =  3)


cangle = seq(0.38,0,length.out=100)
S11 = sin(cangle)*0.2
C11 = cos(cangle)*0.2

cangle = seq(1.03,0,length.out=100)
S12 = sin(cangle)*0.35
C12 = cos(cangle)*0.35

PT9 = PT8+geom_line(aes(x=C11+OBS[1,1],y=S11+OBS[1,2]), linetype = 1,size=0.4)+
geom_line(aes(x=C12+OBS[2,1],y=S12+OBS[2,2]), linetype = 1,size=0.4)+
geom_text(aes(x=0.40,y=0.05), label=deparse(bquote(paste(phi['t'['i']]))),parse=TRUE,size=7)+
geom_text(aes(x=-0.45,y=-0.25), label=deparse(bquote(paste(phi['t'['i-1']]))),parse=TRUE,size=7)
#+
#geom_text(aes(x=-0.1,y=0.9), label=deparse(bquote(paste('v'['t'['i']*','][2]))),parse=TRUE,size=7)+
#geom_text(aes(x=0.5,y=-0.1), label=deparse(bquote(paste('v'['t'['i']*','][1]))),parse=TRUE,size=7)
PT9


library(car)
library(mvtnorm)

app = ellipse(c(0,0), matrix(c(1.9,-0.7,-0.7,0.5), ncol=2)*0.02, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
app[,1] = app[,1]+0.35
app[,2] = app[,2]+1

PT10 = PT9+geom_path(aes(x=app[,1],y=app[,2]),  size = 1)+
geom_segment(aes(x = 0, y = 0, xend = 0.35, yend = 1, color= "Prev. Dir."),arrow = arrow(length = unit(0.3, "cm"), type = "closed"), size = 1, linetype=1, color= "black")+scale_x_continuous(breaks= c(-1,-0.5,0, 0.5,1), labels = c(5,6,7,8,9), limits=c(-1,1))+scale_y_continuous(breaks= c(-1,-0.5,0, 0.5,1), labels = c(5,6,7,8,9)+5, limits=c(-0.7,1.3))
PT10+ coord_fixed()
pdf(paste(DIR_OUT,"TransX.pdf",sep=""),width=7, height=7)
PT10+ coord_fixed()
dev.off()



#### #### #### #### #### #### #### #### ####
####
#### #### #### #### #### #### #### #### ####


Mu = matrix(c(0,0), ncol=1)
Nu = matrix(c(0,6), ncol=1)
tau = 0.25
sigma21 = 0.2
sigma22 = 1
corr   = 0
rho = c(0.33,0.66)
Sigma = matrix(c(sigma21,(sigma21*sigma22)^0.5*corr,(sigma21*sigma22)^0.5*corr,sigma22), ncol=2)

angle = rev(rep(seq(0,-pi,by=-pi/3), each=4))
xseq = yseq = seq(-20,20, length.out=4)
xseq = yseq = seq(-20,20, length.out=4)
Grid = expand.grid(xseq, yseq)
colnames(Grid) = c("Longitude", "Latitude")

# BRW




P = qplot(c(210),0,geom="line", xlab="",ylab="")+theme(axis.text=element_text(size=25),
        axis.title=element_text(size=28,face="bold"))+ coord_fixed(ratio=1)



meanPlot = Grid[,]
meanPlot[,1] = Grid[,1]+tau*(Mu[1,1]-Grid[,1])
meanPlot[,2] = Grid[,2]+tau*(Mu[2,1]-Grid[,2])

P1 = P+geom_point( aes(Grid$Longitude, Grid$Latitude))+geom_point()
#geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = meanPlot[,1], yend = meanPlot[,2]),
                 # arrow = arrow(length = unit(0.2, "cm")),color=cbPalette[1])

meanPlot_2 = meanPlot

i = 1
SigmaPlot = Sigma
app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
app[,1] = app[,1]+meanPlot[i,1]
app[,2] = app[,2]+meanPlot[i,2]
ell = app
for(i in 2:nrow(Grid))
{

    app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
    app[,1] = app[,1]+meanPlot[i,1]
    app[,2] = app[,2]+meanPlot[i,2]
    ell = rbind(ell,NA, app)

}
ell_BRW = ell
P_BWR = P1

# CRW


#P = qplot(c(0),0,geom="line", xlab="",ylab="")+theme(axis.text=element_text(size=25),
        #axis.title=element_text(size=28,face="bold"))+scale_colour_continuous(guide = FALSE)+theme(legend.position="none")+ coord_fixed(ratio=1)

P1 = P_BWR


i = 1
arr = Grid

ang = angle[i]
R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
muPlot = R%*%Nu
arr[i,] = arr[i,]+muPlot
SigmaPlot = R%*%Sigma%*%t(R)
app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
app[,2] = app[,2]+Grid[i,2]+muPlot[2,1]
ell = app
for(i in 2:nrow(Grid))
{
    ang = angle[i]
    R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
    muPlot = R%*%Nu
    SigmaPlot = R%*%Sigma%*%t(R)
    arr[i,] = arr[i,]+muPlot
    app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
    app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
    app[,2] = app[,2]+Grid[i,2]+muPlot[2,1]
    ell = rbind(ell,NA, app)

}

cc = factor(c(1,2,3))



col_ell_1 = rep(cbPalette[1], nrow(ell_BRW))
col_ell_2 = rep(cbPalette[4], nrow(ell))
ell = rbind(ell,NA, ell_BRW)

P4 = P1# + geom_path(aes(x=ell[,1],y=ell[,2], color=c(col_ell_2,NA,col_ell_1)), color=c(col_ell_2,NA,col_ell_1), size = 1)


i = 1
ang = angle[i]
R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
arr2 = Grid
muPlot = Grid[i,]+(R)%*%matrix(c(-8,0), ncol=1)
arr2[i,] = muPlot
for(i in 2:nrow(Grid))
{
    ang = angle[i]
    R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
    muPlot = Grid[i,]+(R)%*%matrix(c(-8,0), ncol=1)
    arr2[i,] = muPlot
}
P6 = P4 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = arr[,1], yend = arr[,2], color="CRW"),
                  arrow = arrow(length = unit(0.3, "cm")), size = 1.2)
P7 = P6 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = meanPlot_2[,1], yend = meanPlot_2[,2], color="BRW"),
                arrow = arrow(length = unit(0.3, "cm")), size = 1.2)

P8 = P7+geom_segment(aes(x = arr2[,1], y = arr2[,2], xend = Grid[,1], yend = Grid[,2], color= "Prev. Dir."),arrow = arrow(length = unit(0.3, "cm"), type = "closed"), size = 1.2, linetype=6, color= "black")

PTOT = P8+scale_color_manual(values = c(cbPalette[c(1,4)],"black"),name="", labels=c(expression(rho ~"=0"),expression(rho ~"=1")))+  theme(axis.text.y = element_text(face="bold",size=25),
axis.text.x = element_text(face="bold",size=25),
axis.title.x = element_text(face="bold",size=25),
axis.title.y = element_text(face="bold",size=25),
legend.text = element_text(face="bold",size=25),
legend.title = element_text(face="bold",size=25) )+
geom_path(aes(x=ell[1:length(col_ell_2),1],y=ell[1:length(col_ell_2),2]),  size = 1, color=c(col_ell_2))+
geom_path(aes(x=ell[-c(1:(length(col_ell_2)+1)),1],y=ell[-c(1:(length(col_ell_2)+1)),2]),  size = 1., color=c(col_ell_1))+
theme(legend.position="bottom")
PTOT+xlim(c(-30,30))+ylim(c(-30,30))

# PTOT =P8+scale_color_manual(values = c(cbPalette[c(1,4)],"black"),name="", labels=c(expression(rho ~"=1/3"),expression(rho ~"=2/3")))+  theme(axis.text.y = element_text(face="bold",size=25),


pdf(paste(DIR_OUT,"EsMov.pdf",sep=""),width=7, height=7)
PTOT+xlim(c(-30,30))+ylim(c(-30,30))
dev.off()




### ### ### ### ###
### STAP
### ### ### ### ###



P = qplot(c(120),0,geom="line", xlab="",ylab="")+theme(axis.text=element_text(size=25),
        axis.title=element_text(size=28,face="bold"))+ coord_fixed(ratio=1)



meanPlot = Grid[,]
meanPlot[,1] = Grid[,1]+tau*(Mu[1,1]-Grid[,1])
meanPlot[,2] = Grid[,2]+tau*(Mu[2,1]-Grid[,2])

P1 = P+geom_point( aes(Grid$Longitude, Grid$Latitude))+geom_point()
#geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = meanPlot[,1], yend = meanPlot[,2]),
                 # arrow = arrow(length = unit(0.2, "cm")),color=cbPalette[1])

j=1
i = 1
arr = Grid
ang = angle[i]
R = matrix(c(cos(rho[j]*ang), sin(rho[j]*ang), -sin(rho[j]*ang), cos(rho[j]*ang)), ncol=2)
R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
muPlot =  (1-rho[j])*tau*(0-Grid[i,])+ rho[j]*R2%*%Nu
arr[i,] = arr[i,]+muPlot
SigmaPlot = R%*%Sigma%*%t(R)
app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
ell = app
for(i in 2:nrow(Grid))
{
    ang = angle[i]
    R = matrix(c(cos(rho[j]*ang), sin(rho[j]*ang), -sin(rho[j]*ang), cos(rho[j]*ang)), ncol=2)
    R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
    muPlot =  (1-rho[j])*tau*(0-Grid[i,])+ rho[j]*R2%*%Nu
    SigmaPlot = R%*%Sigma%*%t(R)
    arr[i,] = arr[i,]+muPlot
    app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
    app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
    app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
    ell = rbind(ell,NA, app)

}
ell_2 = ell
arr_2 = arr
muPlot_2 = arr
j=2
i = 1
arr = Grid
ang = angle[i]
R = matrix(c(cos(rho[j]*ang), sin(rho[j]*ang), -sin(rho[j]*ang), cos(rho[j]*ang)), ncol=2)
R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
muPlot =  (1-rho[j])*tau*(0-Grid[i,])+ rho[j]*R2%*%Nu
arr[i,] = arr[i,]+muPlot
SigmaPlot = R%*%Sigma%*%t(R)
app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
ell = app
for(i in 2:nrow(Grid))
{
    ang = angle[i]
    R = matrix(c(cos(rho[j]*ang), sin(rho[j]*ang), -sin(rho[j]*ang), cos(rho[j]*ang)), ncol=2)
    R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
    muPlot =  (1-rho[j])*tau*(0-Grid[i,])+ rho[j]*R2%*%Nu
    SigmaPlot = R%*%Sigma%*%t(R)
    arr[i,] = arr[i,]+muPlot
    app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
    app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
    app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
    ell = rbind(ell,NA, app)

}

col_ell_1 = rep(cbPalette[1], nrow(ell_2))
col_ell_2 = rep(cbPalette[4], nrow(ell))
ell = rbind(ell,NA, ell_2)

P4 = P1# + geom_path(aes(x=ell[,1],y=ell[,2], color=c(col_ell_2,NA,col_ell_1)), color=c(col_ell_2,NA,col_ell_1), size = 1)


i = 1
ang = angle[i]
R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
arr2 = Grid
muPlot = Grid[i,]+(R)%*%matrix(c(-8,0), ncol=1)
arr2[i,] = muPlot
for(i in 2:nrow(Grid))
{
    ang = angle[i]
    R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
    muPlot = Grid[i,]+(R)%*%matrix(c(-8,0), ncol=1)
    arr2[i,] = muPlot
}
P6 = P4 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = arr[,1], yend = arr[,2], color="2"),
                  arrow = arrow(length = unit(0.3, "cm")), size = 1.2)
P7 = P6 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = muPlot_2[,1], yend = muPlot_2[,2], color="1"),
                arrow = arrow(length = unit(0.3, "cm")), size = 1.2)

P8 = P7+geom_segment(aes(x = arr2[,1], y = arr2[,2], xend = Grid[,1], yend = Grid[,2]),arrow = arrow(length = unit(0.3, "cm"), type = "closed"), size = 1.2, linetype=6, color= "black")

PTOT =P8+scale_color_manual(values = c(cbPalette[c(1,4)],"black"),name="", labels=c(expression(rho ~"=1/3"),expression(rho ~"=2/3")))+  theme(axis.text.y = element_text(face="bold",size=25),
axis.text.x = element_text(face="bold",size=25),
axis.title.x = element_text(face="bold",size=25),
axis.title.y = element_text(face="bold",size=25),
legend.text = element_text(face="bold",size=25),
legend.title = element_text(face="bold",size=25),legend.text.align = 0 )+
geom_path(aes(x=ell[1:length(col_ell_2),1],y=ell[1:length(col_ell_2),2]),  size = 1, color=c(col_ell_2))+
geom_path(aes(x=ell[-c(1:(length(col_ell_2)+1)),1],y=ell[-c(1:(length(col_ell_2)+1)),2]),  size = 1, color=c(col_ell_1))+
theme(legend.position="bottom")
PTOT





pdf(paste(DIR_OUT,"EsMov2.pdf",sep=""),width=7, height=7)
PTOT+xlim(c(-30,30))+ylim(c(-30,30))
dev.off()



# ### ### ### ### ### ### ### ### ###
# ### REAL DATA
# ### ### ### ### ### ### ### ### ###


library( pdfCluster)
#install.packages("pdfCluster")

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
DIR_OUT  = "/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/lavori/HDP_Animal/private_repo/tex/"
doMAPBALLS = F
nanim           = length(MCMCMODEL_OUT)
vecOrdK         = list()
IndexMaxMCMC    = list()
Kmcmc           = list()
ZMAP            = list()
MAP_REV         = list()


table(apply(MCMCMODEL_OUT[[1]]$zeta,1, function(x) length(unique(x))))/2500
table(apply(MCMCMODEL_OUT[[2]]$zeta,1, function(x) length(unique(x))))/2500
table(apply(MCMCMODEL_OUT[[3]]$zeta,1, function(x) length(unique(x))))/2500
table(apply(MCMCMODEL_OUT[[4]]$zeta,1, function(x) length(unique(x))))/2500
table(apply(MCMCMODEL_OUT[[5]]$zeta,1, function(x) length(unique(x))))/2500
table(apply(MCMCMODEL_OUT[[6]]$zeta,1, function(x) length(unique(x))))/2500
### FIgures - Posterior distribution of the number of behaviors - Not in the paper
for(ian in 1:nanim)
{
    print(ian)
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


  if(doMAPBALLS)
  {
      aaaa = comp.psm(ModelOUT$zeta[IndexMaxMCMC[[ian]],])
      app = as.numeric(factor(ZMAP[[ian]], levels=unique(ZMAP[[ian]])[order(unique(ZMAP[[ian]]))]))
      MAP_REV[[ian]] = greedy(aaaa, start.cl = app,suppress.comment=FALSE, L=1000, maxiter=100)
      table(ZMAP[[ian]],MAP_REV[[ian]]$cl)
      save(ZMAP, file = paste(DIR_OUT,ian,"Zmap.Rdata", sep=""))
  }


  TT = table(ZMAP[[ian]])
  vecOrdK[[ian]] =  as.numeric(names(TT[order(-TT)]))
}

if(doMAPBALLS)
{
    for(ian in 1:nanim)
    {
        MAP_REV_app =    ZMAP[[ian]][order(unique(ZMAP[[ian]]))][MAP_REV$cl[[ian]]]
        ZMAP[[ian]] = MAP_REV[[ian]]
        TT = table(ZMAP[[ian]])
        vecOrdK[[ian]] =  as.numeric(names(TT[order(-TT)]))
    }
    save(ZMAP, paste(DIR_OUT,"Zmap.Rdata", sep=""))
}

table(apply(MCMCMODEL_OUT[[1]]$zeta, 1, function(x) length(unique(x))))/2500
table(apply(MCMCMODEL_OUT[[2]]$zeta, 1, function(x) length(unique(x))))/2500
table(apply(MCMCMODEL_OUT[[3]]$zeta, 1, function(x) length(unique(x))))/2500
table(apply(MCMCMODEL_OUT[[4]]$zeta, 1, function(x) length(unique(x))))/2500
table(apply(MCMCMODEL_OUT[[5]]$zeta, 1, function(x) length(unique(x))))/2500
table(apply(MCMCMODEL_OUT[[6]]$zeta, 1, function(x) length(unique(x))))/2500
table(ZMAP[[3]] )/length(ZMAP[[3]] )
#vecOrdK[[3]] = vecOrdK[[3]][c(2,3,4,1)]

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
    dd = paste("$backslash mu_{j,k,",i,"}$ ",sep="")
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
    dd = paste("$backslash eta_{j,k,",i,"}$ ",sep="")
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
    dd = paste("$backslash tau_{j,k}$ ",sep="")
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
    dd = paste("$backslash rho_{j,k}$ ",sep="")
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
      dd = paste("$backslash boldsymbol{backslash Sigma}_{j,k,", i,",",j, "}$ ",sep="")
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
    dd = paste("${backslash pi}_{j,", ik,",","k", "}$ ",sep="")
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
  dd = paste("$n_{j,k}$ ",sep="")
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

table(ZMAP[[3]])

RANDINDEX = matrix(NA, ncol=6, nrow=6)
for(i in 1:6)
{
   for(j in 1:6)
   {
      RANDINDEX[i,j] = adj.rand.index(ZMAP[[i]],ZMAP[[j]])
   }
}
plot(c(RANDINDEX))

round(RANDINDEX,3)
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

SimMatMu = matrix(0, ncol=sum(unlist(Kmcmc)),nrow=sum(unlist(Kmcmc)))
SimMatEta = matrix(0, ncol=sum(unlist(Kmcmc)),nrow=sum(unlist(Kmcmc)))
SimMatNu = matrix(0, ncol=sum(unlist(Kmcmc)),nrow=sum(unlist(Kmcmc)))
SimMatRho = matrix(0, ncol=sum(unlist(Kmcmc)),nrow=sum(unlist(Kmcmc)))
SimMatSigma = matrix(0, ncol=sum(unlist(Kmcmc)),nrow=sum(unlist(Kmcmc)))



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

diag(SimMatEta) = 1
diag(SimMatSigma) = 1
diag(SimMatMu) = 1
diag(SimMatSigma) = 1
diag(SimMatNu) = 1
diag(SimMatRho) = 1
# # Index = c(
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
# i=1
# for(i in 1:5)
# {
#   if(i==1){ddd = SimMatMu}
#   if(i==2){ddd = SimMatEta}
#   if(i==3){ddd = SimMatRho}
#   if(i==4){ddd = SimMatNu}
#   if(i==5){ddd = SimMatSigma}
#   dataGG = melt(ddd)
#   dataGG$Var2 = max(dataGG$Var2)-dataGG$Var2+1
#
#   P = ggplot(dataGG , aes(Var1,Var2, fill=value)) + geom_raster()+scale_fill_gradient2(low="white",mid="#b2e2e2", high="#2c7fb8",limits=c(0,1), midpoint = (0.5))+
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank(),
#     #axis.line = element_line(colour = "black"),
#     axis.title.x=element_blank(),
#     axis.text.x=element_blank(),
#     legend.title = element_text( size = 24),
#     legend.text = element_text( size = 20),
#     axis.ticks.x=element_blank(),
#     axis.title.y=element_blank(),
#     axis.text.y=element_blank(),
#     axis.ticks.y=element_blank()
#   )+labs(fill = "Prob.")+
#   geom_hline(yintercept=1:(nrow(SimMatNu)+1)-0.5,size=0.1)+
#   geom_vline(xintercept=(1:(nrow(SimMatNu)+1)-0.5) ,size=0.1)+
#   geom_hline(yintercept=nrow(SimMatNu)+1-(c(1,cumsum(unlist(Kmcmc))+1)-0.5),size=1.1)+
#   geom_vline(xintercept=((c(1,cumsum(unlist(Kmcmc))+1)-0.5)),size=1.1)+
#   coord_fixed(ratio = 1)+annotate("text", x = 1:19, y=20.3, label = c(
#   "B11","B12","B13",
#   "B21","B22","B23",
#   "B31","B32","B33","B34",
#   "B41","B42","B43",
#   "B51","B52","B53",
#   "B61","B62","B63"
#
#   ),angle = -90, size = 5
# )+annotate("text", y = rev(1:19), x=-0.3, label = c(
#   "B11","B12","B13",
#   "B21","B22","B23",
#   "B31","B32","B33","B34",
#   "B41","B42","B43",
#   "B51","B52","B53",
#   "B61","B62","B63"
#
#   ),angle = 0, size = 5
# )+annotate("text", y = rev(c(2,5,8,11.5,15,18)), x=-2, label = c(
#   "Woody", "Sherl.","Alvin","Rosie","Bear","Lucy"
#
#   ),angle = 90, size = 5
# )+annotate("text", x = (c(2,5,8.5,12,15,18)), y=20+2, label = c(
#   "Woody", "Sherl.","Alvin","Rosie","Bear","Lucy"
#
#   ),angle = 0, size = 5
# )
# P
#
#   png(paste(DIR_OUT ,"Sim_",NAMEPAR[i],".png",sep=""))
#   print(P)
#   dev.off()
# }



SimMatEta  = SimMatEta[-13,-13]
SimMatSigma = SimMatSigma[-13,-13]
SimMatMu = SimMatMu[-13,-13]
SimMatNu = SimMatNu [-13,-13]
SimMatRho = SimMatRho[-13,-13]


SimMatEta  = SimMatEta[-10,-10]
SimMatSigma = SimMatSigma[-10,-10]
SimMatMu = SimMatMu[-10,-10]
SimMatNu = SimMatNu[-10,-10]
SimMatRho = SimMatRho[-10,-10]


Kmcmc2 = c(3,3,3,2,3,3)
i=1
for(i in 1:5)
{
  if(i==1){ddd = SimMatMu}
  if(i==2){ddd = SimMatEta}
  if(i==3){ddd = SimMatRho}
  if(i==4){ddd = SimMatNu}
  if(i==5){ddd = SimMatSigma}


#"#b2e2e2"
  #ddd[(ddd)<0.25] = NA
  dataGG = melt(ddd)
  dataGG$Var2 = max(dataGG$Var2)-dataGG$Var2+1

  P = ggplot(dataGG , aes(Var1,Var2, fill=value)) + geom_raster()+scale_fill_gradient2(low="white",mid="#b2e2e2", high="red",limits=c(0,1), midpoint = (0.60),na.value="grey")+
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
  geom_hline(yintercept=nrow(SimMatNu)+1-(c(1,cumsum(unlist(Kmcmc2))+1)-0.5),size=1.1)+
  geom_vline(xintercept=((c(1,cumsum(unlist(Kmcmc2))+1)-0.5)),size=1.1)+
  coord_fixed(ratio = 1)+annotate("text", x = 1:17, y=20.3-2, label = c(
      "B[11]","B[12]","B[13]",
      "B[21]","B[22]","B[23]",
      "B[31]","B[32]","B[33]",
      "B[41]","B[42]",
      "B[51]","B[52]","B[53]",
      "B[61]","B[62]","B[63]"

  ),angle = -90, size = 5, parse=TRUE
)+annotate("text", y = rev(1:17), x=-0.3, label = c(
  "B[11]","B[12]","B[13]",
  "B[21]","B[22]","B[23]",
  "B[31]","B[32]","B[33]",
  "B[41]","B[42]",
  "B[51]","B[52]","B[53]",
  "B[61]","B[62]","B[63]"

  ),angle = 0, size = 5, parse=TRUE
)+annotate("text", y = rev(c(2,5,8-0.5,11.5-1.5,15-2,18-2)), x=-1.5, label = c(
  "Woody", "Sherl.","Alvin","Rosie","Bear","Lucy"

),angle = 90, size = 4.9
)+annotate("text", x = (c(2,5,8,12-1.5,15-2,18-2)), y=20+2-2.5, label = c(
  "Woody", "Sherl.","Alvin","Rosie","Bear","Lucy"

),angle = 0, size = 4.9
)
P

  png(paste(DIR_OUT ,"Sim_",NAMEPAR[i],".png",sep=""))
  print(P)
  dev.off()
}


#
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
#     legend.title = element_text( size = 24),
#     legend.text = element_text( size = 20),
#     axis.text.x=element_blank(),
#     axis.ticks.x=element_blank(),
#     axis.title.y=element_blank(),
#     axis.text.y=element_blank(),
#     axis.ticks.y=element_blank()
#   )+labs(fill = "Prob.")+
#   geom_hline(yintercept=1:(nrow(SimMatNu)+1)-0.5,size=0.1)+
#   geom_vline(xintercept=(1:(nrow(SimMatNu)+1)-0.5) ,size=0.1)+
#   geom_hline(yintercept=nrow(SimMatNu)+1-(c(1,cumsum(unlist(Kmcmc))+1)-0.5),size=1.1)+
#   geom_vline(xintercept=((c(1,cumsum(unlist(Kmcmc))+1)-0.5)),size=1.1)+
#   coord_fixed(ratio = 1)+annotate("text", x = 1:19, y=20.3, label = c(
#   "B11","B12","B13",
#   "B21","B22","B23",
#   "B31","B32","B33","B34",
#   "B41","B42","B43",
#   "B51","B52","B53",
#   "B61","B62","B63"
#
#   ),angle = -90, size = 5
# )+annotate("text", y = rev(1:19), x=-0.3, label = c(
#   "B11","B12","B13",
#   "B21","B22","B23",
#   "B31","B32","B33","B34",
#   "B41","B42","B43",
#   "B51","B52","B53",
#   "B61","B62","B63"
#
#   ),angle = 0, size = 5
# )+annotate("text", y = rev(c(2,5,8,11.5,15,18)), x=-2, label = c(
#   "Woody", "Sherl.","Alvin","Rosie","Bear","Lucy"
#
#   ),angle = 90, size = 5
# )+annotate("text", x = (c(2,5,8.5,12,15,18)), y=20+2, label = c(
#   "Woody", "Sherl.","Alvin","Rosie","Bear","Lucy"
#
#   ),angle = 0, size = 5
# )
#
#
#   png(paste(DIR_OUT ,"Sim_Proc",".png",sep=""))
#   print(P)
#   dev.off()


#### plot behaviours
str(Kmcmc)
ianim = 3
for(ianim in 1:nanim)
{

    mu0List = list()
    muCList = list()
    psiList = list()
    rhoList = list()
    sigmaList = list()
    piList = list()
    kk = 1
    ModelOUT = MCMCMODEL_OUT[[ianim]]
    WM         =  vecOrdK[[ianim]][1:unlist(Kmcmc)[ianim]]
    for(k in WM)
    {
        mu0List[[kk]] = matrix(colMeans(ModelOUT$mu[,1:2,k]),ncol=1)
        muCList[[kk]] = matrix(colMeans(ModelOUT$eta[,1:2,k]),ncol=1)
        psiList[[kk]] = matrix(mean(ModelOUT$nu[,1,k]),ncol=1)
        rhoList[[kk]] = matrix(mean(ModelOUT$rho[,1,k]),ncol=1)
        sigmaList[[kk]] = matrix(colMeans(ModelOUT$sigma[,1:4,k]),nrow=2)
        piList[[kk]] = matrix(colMeans(ModelOUT$pi[,1:length(WM),k]),ncol=length(WM))
        kk = kk+1

    }
    kk = 1
    for(k in WM)
    {
        mu0List[[kk]] = matrix(ModelOUT$mu[IndexMAP [ianim],1:2,k],ncol=1)
        muCList[[kk]] = matrix(ModelOUT$eta[IndexMAP [ianim],1:2,k],ncol=1)
        psiList[[kk]] = matrix(ModelOUT$nu[IndexMAP [ianim],1,k],ncol=1)
        rhoList[[kk]] = matrix(ModelOUT$rho[IndexMAP [ianim],1,k],ncol=1)
        sigmaList[[kk]] = matrix(ModelOUT$sigma[IndexMAP [ianim],1:4,k],nrow=2)
        piList[[kk]] = matrix(ModelOUT$pi[IndexMAP [ianim],1:length(WM),k],ncol=length(WM))
        kk = kk+1
    }
    kk = 1
    for(k in WM)
    {
        mu0List[[kk]] = matrix(ModelOUT$mu[IndexMaxDensList[[ianim]][k],1:2,k],ncol=1)
        muCList[[kk]] = matrix(ModelOUT$eta[IndexMaxDensList[[ianim]][k],1:2,k],ncol=1)
        psiList[[kk]] = matrix(ModelOUT$nu[IndexMaxDensList[[ianim]][k],1,k],ncol=1)
        rhoList[[kk]] = matrix(ModelOUT$rho[IndexMaxDensList[[ianim]][k],1,k],ncol=1)
        sigmaList[[kk]] = matrix(ModelOUT$sigma[IndexMaxDensList[[ianim]][k],1:4,k],nrow=2)
        piList[[kk]] = matrix(ModelOUT$pi[IndexMaxDensList[[ianim]][k],1:length(WM),k],ncol=length(WM))
        kk = kk+1

    }
    # for(k in WM)
    # {
    #     mu0List[[kk]] = matrix(apply(ModelOUT$mu[,1:2,k],2,median),ncol=1)
    #     muCList[[kk]] = matrix(apply(ModelOUT$eta[,1:2,k],2,median),ncol=1)
    #     psiList[[kk]] = matrix(median(ModelOUT$nu[,1,k]),ncol=1)
    #     rhoList[[kk]] = matrix(median(ModelOUT$rho[,1,k]),ncol=1)
    #     sigmaList[[kk]] = matrix(apply(ModelOUT$sigma[,1:4,k],2,median),nrow=2)
    #     piList[[kk]] = matrix(apply(ModelOUT$pi[,1:length(WM),k],2,median),ncol=length(WM))
    #     kk = kk+1
    #
    # }

    ### GIAN sistemare da qui
#    DataZ = data.frame(Longitude = DataCoords[,(ianin-1)*2+1], Latitude = DataCoords[,(ianin-1)*2+2] )

    DataZ = data.frame(Longitude = DataCoords[-nrow(DataCoords),(ianim-1)*2+1], Latitude = DataCoords[-nrow(DataCoords),(ianim-1)*2+2],Cluster = as.factor(ZMAP[[ianim]]) )


    kk = 1
    for(kk in 1:length(WM))
    {

        cccol = c(1,2,3,4)[kk]
        xxlim = c(0.09,1,2.7,2.4)[kk]
        xxlim2 = (c(0.09,1,2.7,2.4)*1.6)[kk]
        if((ianim==3))
        {
            xxlim = c(2.7,1,0.09,2.7)[kk]
            xxlim2 = (c(2.7,0.9,0.09,2.7)*1.6)[kk]

            cccol = c(3,2,1,3)[kk]
        }
        Mu = mu0List[[kk]]
        Nu = muCList[[kk]]
        tau = psiList[[kk]]


        rho = rhoList[[kk]]
        Sigma = sigmaList[[kk]]

        xseq = yseq = seq(-xxlim,xxlim,length.out=3)
        xseq = yseq = seq(-xxlim,xxlim, length.out=3)
        nnnn = length(xseq)
        angle = rev(rep(seq(0,-pi,by=-pi/(nnnn-1)), each=nnnn))

        Grid = expand.grid(xseq, yseq)
        colnames(Grid) = c("Longitude", "Latitude")


        P =qplot(c(-10),0,geom="line", xlab="",ylab="")+theme(axis.text.x = element_text(face="bold",size=25),
                axis.text.y = element_text(face="bold",size=25),
                axis.title.x = element_text(face="bold",size=25),
                axis.title.y = element_text(face="bold",size=25),)+ coord_fixed(ratio=1)



        meanPlot = Grid[,]
        meanPlot[,1] = Grid[,1]+tau*(Mu[1,1]-Grid[,1])
        meanPlot[,2] = Grid[,2]+tau*(Mu[2,1]-Grid[,2])

        P1 = P+geom_point( aes(Grid$Longitude, Grid$Latitude))+geom_point()
        #geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = meanPlot[,1], yend = meanPlot[,2]),
                         # arrow = arrow(length = unit(0.2, "cm")),color=cbPalette[1])

        j=1
        i = 1
        arr = Grid
        ang = angle[i]
        R = matrix(c(cos(rho*ang), sin(rho*ang), -sin(rho*ang), cos(rho*ang)), ncol=2)
        R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
        muPlot =  (1-rho[j])*tau*(Mu-Grid[i,])+ rho[j]*R2%*%Nu
        arr[i,] = arr[i,]+muPlot
        SigmaPlot = R%*%Sigma%*%t(R)
        app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
        app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
        app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
        ell = app
        for(i in 2:nrow(Grid))
        {
            ang = angle[i]
            R = matrix(c(cos(rho*ang), sin(rho*ang), -sin(rho*ang), cos(rho*ang)), ncol=2)
            R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
            muPlot =  (1-rho[j])*tau*(Mu-Grid[i,])+ rho[j]*R2%*%Nu
            SigmaPlot = R%*%Sigma%*%t(R)
            arr[i,] = arr[i,]+muPlot
            app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
            app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
            app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
            ell = rbind(ell,NA, app)

        }
        ell_2 = ell
        arr_2 = arr
        muPlot_2 = arr


        # col_ell_1 = rep(cbPalette[1], nrow(ell_2))
        # col_ell_2 = rep(cbPalette[3], nrow(ell))
        # ell = rbind(ell,NA, ell_2)

        P4 = P1



        i = 1
        ang = angle[i]
        R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
        arr2 = Grid
        muPlot = Grid[i,]+(R)%*%matrix(c(-xxlim2/5,0), ncol=1)
        arr2[i,] = muPlot
        for(i in 2:nrow(Grid))
        {
            ang = angle[i]
            R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
            muPlot = Grid[i,]+(R)%*%matrix(c(-xxlim2/5,0), ncol=1)
            arr2[i,] = muPlot
        }
        #P6 = P4 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = arr[,1], yend = arr[,2], color="2"),
                          #arrow = arrow(length = unit(0.2, "cm")))
        #P7 = P6 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = muPlot_2[,1], yend = muPlot_2[,2], color="1"),
                        #arrow = arrow(length = unit(0.2, "cm")))


        P5 = P4+geom_segment(aes(x = arr2[,1], y = arr2[,2], xend = Grid[,1], yend = Grid[,2]),arrow = arrow(length = unit(0.3, "cm"), type = "closed"), size = 1.2, linetype=6, color= "black")
        P7 = P5 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = muPlot_2[,1], yend = muPlot_2[,2]), color=cbPalette[c(cccol)],
                        arrow = arrow(length = unit(0.3, "cm")), size = 1.2)

        PTOT =P7+
        #scale_color_manual(values = c(cbPalette[c(kk)],"black",cbPalette[c(kk)]),name="", labels=c("Behavior","Prev. Dir"))+
         theme(
             axis.text.y = element_text(face="bold",size=25),
        axis.text.x = element_text(face="bold",size=25),
        axis.title.x = element_text(face="bold",size=25),
        axis.title.y = element_text(face="bold",size=25),
        legend.text = element_text(face="bold",size=25),
        legend.title = element_text(face="bold",size=25) )+ geom_path(aes(x=ell_2[,1],y=ell_2[,2]),  size = 2,color=cbPalette[c(cccol)], linetype=1)+ylab("Latitude")+xlab("Longitude")+theme(legend.position="bottom")

        PP = PTOT+xlim(c(-xxlim2,xxlim2))+ylim(c(-xxlim2,xxlim2))
        #+geom_point(aes(DataZ$Longitude[as.numeric(as.character(DataZ$Cluster))==WM[kk]],DataZ$Latitude[as.numeric(as.character(DataZ$Cluster))==WM[kk]]),size=0.05)


    	PP2 = PP

        if((rho<0.9)& (tau>0.1))
        {
            PP2 = PP2+geom_point(aes(Mu[1],Mu[2]),shape = 8, colour = "brown", fill = "white", size = 12, stroke = 3)
        }
    	# if((kk==4) | (kk==5))
    	# {
    	#   Xp = mean(ModelOUT$mu0[,1,WM[i]])
    	#   Yp = mean(ModelOUT$mu0[,2,WM[i]])
    	#   PP2 = PP2 +annotate("point", x = Xp, y = Yp, colour = "black", size=5)
    	#   #+
    	#   #annotate("text", x = Xp, y = Yp+0.4, label="Attractive-point", size=15)
    	# }
    	#PP2

        # if(rho<0.4)
        # {
        #     PP2 = PP2+geom_point(aes(Mu[1],Mu[2]),shape = 8, colour = "black", fill = "white", size = 10, stroke = 3)
        # }

        pdf(paste(DIR_OUT,ianim,"_",kk,"Mov.pdf",sep=""),width=7, height=7)
        print(PP2+xlim(c(-xxlim2,xxlim2))+ylim(c(-xxlim2,xxlim2)))
        dev.off()

    }

}


# ian1 = 1
# ian2 = 6
# table(ZMAP[[ian1]],ZMAP[[ian2]])
# vecOrdK[[ian1]][1:unlist(Kmcmc)[ian1]]
# vecOrdK[[ian2]][1:unlist(Kmcmc)[ian2]]








#### #### #### #### #### ####
#### POSTERIOR PREDICTIVE DISTRIBUTION
#### #### #### #### #### ####

j = 1
str(Kmcmc)
ianim = 3
library(rmutil)
library(matrixcalc)
library(GenFunc)
library(emdbook)
library(coda)
#install.packages("emdbook")

for(ianim in 1:nanim)
{

    mu0List = list()
    muCList = list()
    psiList = list()
    rhoList = list()
    sigmaList = list()
    piList = list()
    kk = 1
    ModelOUT = MCMCMODEL_OUT[[ianim]]
    WM         =  vecOrdK[[ianim]][1:unlist(Kmcmc)[ianim]]

    nsim = dim(ModelOUT$mu)[1]

	 kk = 1
	 cccol = c(1,2,3,4)[kk]
	 xxlim = c(0.09,1,2.7,2.4)[kk]
	 xxlim2 = (c(0.09,1,2.7,2.4)*1.6)[kk]
	 if((ianim==3))
	 {
		  xxlim = c(2.7,1,0.09,2.7)[kk]
		  xxlim2 = (c(2.7,0.9,0.09,2.7)*1.6)[kk]

		  cccol = c(3,2,1,3)[kk]
	 }


	 xseq = yseq = seq(-xxlim,xxlim,length.out=3)
	 xseq = yseq = seq(-xxlim,xxlim, length.out=3)
	 nnnn = length(xseq)
	 angle = rev(rep(seq(0,-pi,by=-pi/(nnnn-1)), each=nnnn))

	 np = length(xseq)^2

	 results = list()
	 for(ii in 1:np)
	 {
		 results[[ii]] = matrix(NA, ncol=2, nrow=nsim)
	 }
	 isim = 1
    for(isim in 1:nsim)
    {
		 kk = 1
		 for(k in WM)
		 {
			  mu0List[[kk]] = matrix(ModelOUT$mu[isim,1:2,k],ncol=1)
			  muCList[[kk]] = matrix(ModelOUT$eta[isim,1:2,k],ncol=1)
			  psiList[[kk]] = matrix(ModelOUT$nu[isim,1,k],ncol=1)
			  rhoList[[kk]] = matrix(ModelOUT$rho[isim,1,k],ncol=1)
			  sigmaList[[kk]] = matrix(ModelOUT$sigma[isim,1:4,k],nrow=2)
			  piList[[kk]] = matrix(ModelOUT$pi[isim,WM,k],ncol=length(WM))
			  kk = kk+1
		 }

		 ## stationary distribution
		 P = matrix(unlist(piList), ncol=length(WM), byrow=T)


		 r=eigen(P)
		 rvec=r$vectors
		 lvec=ginv(r$vectors)
		 lam<-r$values
		 rvec%*%diag(lam)%*%ginv(rvec)
		 Wmin = which.min(abs(lam-1))
		 pi_eig<-lvec[Wmin,]/sum(lvec[Wmin,])


		 kk = sample(1:length(WM),1,prob=pi_eig)


		 Mu = mu0List[[kk]]
		 Nu = muCList[[kk]]
		 tau = psiList[[kk]]


		 rho = rhoList[[kk]]
		 Sigma = sigmaList[[kk]]

		 Grid = expand.grid(xseq, yseq)
		 colnames(Grid) = c("Longitude", "Latitude")

		 i = 1
		 j = 1
		 for(i in 1:np)
		 {
			 ang = angle[i]
			 R = matrix(c(cos(rho*ang), sin(rho*ang), -sin(rho*ang), cos(rho*ang)), ncol=2)
			 R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
			 muPlot =  Grid[i,]+(1-rho[j])*tau*(Mu-Grid[i,])+ rho[j]*R2%*%Nu
			 SigmaPlot = R%*%Sigma%*%t(R)

			 results[[i]][isim,] = c(rmnorm(1,muPlot,SigmaPlot))
		 }

	 }
	 plot(results[[1]], col=1, xlim=c(-3,3), ylim=c(-3,3))
	 points(results[[2]], col=2, xlim=c(-3,3), ylim=c(-3,3))
	 points(results[[3]], col=3, xlim=c(-3,3), ylim=c(-3,3))
	 points(results[[4]], col=4, xlim=c(-3,3), ylim=c(-3,3))
	 points(results[[5]], col=1, xlim=c(-3,3), ylim=c(-3,3))
	 points(results[[6]], col=1, xlim=c(-3,3), ylim=c(-3,3))
	 points(results[[7]], col=1, xlim=c(-3,3), ylim=c(-3,3))
	 points(results[[8]], col=1, xlim=c(-3,3), ylim=c(-3,3))
	 points(results[[9]], col=1, xlim=c(-3,3), ylim=c(-3,3))
	 a = HPDregionplot(as.mcmc(results[[9]][1:2]), n=100, h=0.2, add=T, col=c(2,2,2), lwd=3)

	 plot(results[[1]])
	 Grid
	 a = HPDregionplot(results[[1]][,1:2], add=F, col=c(2,2,2), lwd=3)




    DataZ = data.frame(Longitude = DataCoords[-nrow(DataCoords),(ianim-1)*2+1], Latitude = DataCoords[-nrow(DataCoords),(ianim-1)*2+2],Cluster = as.factor(ZMAP[[ianim]]) )


    kk = 1
    for(kk in 1:length(WM))
    {






        P =qplot(c(-10),0,geom="line", xlab="",ylab="")+theme(axis.text.x = element_text(face="bold",size=25),
                axis.text.y = element_text(face="bold",size=25),
                axis.title.x = element_text(face="bold",size=25),
                axis.title.y = element_text(face="bold",size=25),)+ coord_fixed(ratio=1)



        meanPlot = Grid[,]
        meanPlot[,1] = Grid[,1]+tau*(Mu[1,1]-Grid[,1])
        meanPlot[,2] = Grid[,2]+tau*(Mu[2,1]-Grid[,2])

        P1 = P+geom_point( aes(Grid$Longitude, Grid$Latitude))+geom_point()
        #geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = meanPlot[,1], yend = meanPlot[,2]),
                         # arrow = arrow(length = unit(0.2, "cm")),color=cbPalette[1])

        j=1
        i = 1
        arr = Grid
        ang = angle[i]
        R = matrix(c(cos(rho*ang), sin(rho*ang), -sin(rho*ang), cos(rho*ang)), ncol=2)
        R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
        muPlot =  (1-rho[j])*tau*(Mu-Grid[i,])+ rho[j]*R2%*%Nu
        arr[i,] = arr[i,]+muPlot
        SigmaPlot = R%*%Sigma%*%t(R)
        app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
        app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
        app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
        ell = app
        for(i in 2:nrow(Grid))
        {
            ang = angle[i]
            R = matrix(c(cos(rho*ang), sin(rho*ang), -sin(rho*ang), cos(rho*ang)), ncol=2)
            R2 = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
            muPlot =  (1-rho[j])*tau*(Mu-Grid[i,])+ rho[j]*R2%*%Nu
            SigmaPlot = R%*%Sigma%*%t(R)
            arr[i,] = arr[i,]+muPlot
            app = ellipse(c(0,0), SigmaPlot, col="red", radius=sqrt(2 * qf(.95, 2, 9999)), add=F, draw=F)
            app[,1] = app[,1]+Grid[i,1]+muPlot[1,1]
            app[,2] = app[,2]+Grid[i,2]+muPlot[1,2]
            ell = rbind(ell,NA, app)

        }
        ell_2 = ell
        arr_2 = arr
        muPlot_2 = arr


        # col_ell_1 = rep(cbPalette[1], nrow(ell_2))
        # col_ell_2 = rep(cbPalette[3], nrow(ell))
        # ell = rbind(ell,NA, ell_2)

        P4 = P1



        i = 1
        ang = angle[i]
        R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
        arr2 = Grid
        muPlot = Grid[i,]+(R)%*%matrix(c(-xxlim2/5,0), ncol=1)
        arr2[i,] = muPlot
        for(i in 2:nrow(Grid))
        {
            ang = angle[i]
            R = matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol=2)
            muPlot = Grid[i,]+(R)%*%matrix(c(-xxlim2/5,0), ncol=1)
            arr2[i,] = muPlot
        }
        #P6 = P4 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = arr[,1], yend = arr[,2], color="2"),
                          #arrow = arrow(length = unit(0.2, "cm")))
        #P7 = P6 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = muPlot_2[,1], yend = muPlot_2[,2], color="1"),
                        #arrow = arrow(length = unit(0.2, "cm")))


        P5 = P4+geom_segment(aes(x = arr2[,1], y = arr2[,2], xend = Grid[,1], yend = Grid[,2]),arrow = arrow(length = unit(0.3, "cm"), type = "closed"), size = 1.2, linetype=6, color= "black")
        P7 = P5 + geom_segment(aes(x = Grid[,1], y = Grid[,2], xend = muPlot_2[,1], yend = muPlot_2[,2]), color=cbPalette[c(cccol)],
                        arrow = arrow(length = unit(0.3, "cm")), size = 1.2)

        PTOT =P7+
        #scale_color_manual(values = c(cbPalette[c(kk)],"black",cbPalette[c(kk)]),name="", labels=c("Behavior","Prev. Dir"))+
         theme(
             axis.text.y = element_text(face="bold",size=25),
        axis.text.x = element_text(face="bold",size=25),
        axis.title.x = element_text(face="bold",size=25),
        axis.title.y = element_text(face="bold",size=25),
        legend.text = element_text(face="bold",size=25),
        legend.title = element_text(face="bold",size=25) )+ geom_path(aes(x=ell_2[,1],y=ell_2[,2]),  size = 2,color=cbPalette[c(cccol)], linetype=1)+ylab("Latitude")+xlab("Longitude")+theme(legend.position="bottom")

        PP = PTOT+xlim(c(-xxlim2,xxlim2))+ylim(c(-xxlim2,xxlim2))
        #+geom_point(aes(DataZ$Longitude[as.numeric(as.character(DataZ$Cluster))==WM[kk]],DataZ$Latitude[as.numeric(as.character(DataZ$Cluster))==WM[kk]]),size=0.05)


    	PP2 = PP

        if((rho<0.9)& (tau>0.1))
        {
            PP2 = PP2+geom_point(aes(Mu[1],Mu[2]),shape = 8, colour = "brown", fill = "white", size = 12, stroke = 3)
        }
    	# if((kk==4) | (kk==5))
    	# {
    	#   Xp = mean(ModelOUT$mu0[,1,WM[i]])
    	#   Yp = mean(ModelOUT$mu0[,2,WM[i]])
    	#   PP2 = PP2 +annotate("point", x = Xp, y = Yp, colour = "black", size=5)
    	#   #+
    	#   #annotate("text", x = Xp, y = Yp+0.4, label="Attractive-point", size=15)
    	# }
    	#PP2

        # if(rho<0.4)
        # {
        #     PP2 = PP2+geom_point(aes(Mu[1],Mu[2]),shape = 8, colour = "black", fill = "white", size = 10, stroke = 3)
        # }

        pdf(paste(DIR_OUT,ianim,"_",kk,"Mov.pdf",sep=""),width=7, height=7)
        print(PP2+xlim(c(-xxlim2,xxlim2))+ylim(c(-xxlim2,xxlim2)))
        dev.off()

    }

}
