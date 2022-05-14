# Exam 1

# Q1:
> #E1.r
> 
> rm(list=ls())       # cleans workspace
> getwd()
[1] "/Users/Jing"
> setwd("/Users/Jing/Desktop")
> d0=read.csv("E1.csv",header=T)
> d1=d0[,-c(1:10)]
> head(d1)
   DJIAadj FTSE100adj CAC40adj Nikkeiadj
1 11219.38   11131.84 6373.894  131.7744
2 11173.59   11096.28 6378.162  134.3818
3 11076.18   11185.35 6474.040  135.9433
4 11124.37   11016.71 6357.486  135.4381
5 11088.02   11040.73 6364.764  134.1003
6 11097.87   11109.50 6431.668  136.1710
> n=nrow(d0)
> n
[1] 501
> 
> aux=d1$DJIAadj
> change=diff(aux)     # succesive differences
> row1=aux[-n]       # prices ignoring last row
> Dret=change/row1     # returns
> 
> aux=d1$FTSE100adj
> change=diff(aux)     # succesive differences
> row1=aux[-n]       # prices ignoring last row
> Fret=change/row1     # returns
> 
> aux=d1$CAC40adj
> change=diff(aux)     # succesive differences
> row1=aux[-n]       # prices ignoring last row
> Cret=change/row1     # returns
> 
> aux=d1$Nikkeiadj
> change=diff(aux)     # succesive differences
> row1=aux[-n]       # prices ignoring last row
> Nret=change/row1     # returns
> 
> d2=data.frame(Dret,Fret,Cret,Nret)
> 
> means=colMeans(d2)    # daily mean returns
> means
         Dret          Fret          Cret          Nret 
 2.599304e-05 -1.957257e-04  4.197056e-05 -2.144309e-04 
> 
> cova=cov(d2)          # cov matrix daily returns
> cova
              Dret         Fret         Cret          Nret
Dret  1.229524e-04 7.696591e-05 7.682514e-05 -9.493488e-06
Fret  7.696591e-05 2.013973e-04 1.821076e-04  3.944302e-05
Cret  7.682514e-05 1.821076e-04 1.953506e-04  4.078130e-05
Nret -9.493488e-06 3.944302e-05 4.078130e-05  1.913129e-04
> 
> vars=diag(cova)
> sqrt(vars)
      Dret       Fret       Cret       Nret 
0.01108839 0.01419145 0.01397679 0.01383159


#Volatility:
From Excel, we can calculate the r_g (geo mean) of each index.
DJIAadj		FTSE100adj	  CAC40adj		Nikkeiadj	
0.999964513		0.99970393		0.999944807		0.99968949

s_g=sqrt(1/(n-1)* sum((Dret-0.999964513)^2))
DJIA:
> sqrt(252)*s_g1
FTSE:
> sqrt(252)*s_g2
CAC:
> sqrt(252)*s_g3
Nikkei:
> sqrt(252)*s_g4

# q2
> #efficient portfolio
> #===========================================================
> mu=max(means)
> mu
[1] 4.197056e-05

> x=c(-means)
> x
         Dret          Fret          Cret          Nret 
-2.599304e-05  1.957257e-04 -4.197056e-05  2.144309e-04 
> y=c(-1,-1,-1,-1)
> aux=cbind(cova,x,y)
> aux
              Dret         Fret         Cret          Nret             x  y
Dret  1.229524e-04 7.696591e-05 7.682514e-05 -9.493488e-06 -2.599304e-05 -1
Fret  7.696591e-05 2.013973e-04 1.821076e-04  3.944302e-05  1.957257e-04 -1
Cret  7.682514e-05 1.821076e-04 1.953506e-04  4.078130e-05 -4.197056e-05 -1
Nret -9.493488e-06 3.944302e-05 4.078130e-05  1.913129e-04  2.144309e-04 -1

> s=c(means,0,0)
> t=c(1,1,1,1,0,0)
> a=rbind(aux,s,t)
> a
              Dret          Fret         Cret          Nret             x  y
Dret  1.229524e-04  7.696591e-05 7.682514e-05 -9.493488e-06 -2.599304e-05 -1
Fret  7.696591e-05  2.013973e-04 1.821076e-04  3.944302e-05  1.957257e-04 -1
Cret  7.682514e-05  1.821076e-04 1.953506e-04  4.078130e-05 -4.197056e-05 -1
Nret -9.493488e-06  3.944302e-05 4.078130e-05  1.913129e-04  2.144309e-04 -1
s     2.599304e-05 -1.957257e-04 4.197056e-05 -2.144309e-04  0.000000e+00  0
t     1.000000e+00  1.000000e+00 1.000000e+00  1.000000e+00  0.000000e+00  0

> b=c(0,0,0,0,mu,1)
> b
[1] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 4.197056e-05 1.000000e+00

> weight=solve(a,b)
> weight
         Dret          Fret          Cret          Nret             x             y 
 6.075003e-01 -4.034189e-01  4.597863e-01  3.361323e-01  5.988776e-02  7.421956e-05
> weight=weight[-c(5,6)]
> weight
      Dret       Fret       Cret       Nret 
 0.6075003 -0.4034189  0.4597863  0.3361323 

> meanfound=t(weight)%*%means
> meanfound
             [,1]
[1,] 4.197056e-05

> var=t(weight)%*%cova%*%weight
> var
             [,1]
[1,] 7.673308e-05
> sqrt(var)
            [,1]
[1,] 0.008759742

> R=0.6075003*Dret-0.4034189*Fret+0.4597863*Cret+0.3361323*Nret
r_g=GEOMEAN(1+R)
s_g=sqrt(1/(n-1)*sum((Dret-0.999964513)^2))
annual volatility=sqrt(252)*s_g

