library(wavelets)

#Threshold function
#
#  Function k=Threshold(alpha,n)
#
#  This function computes the threshold value to be used in the programs
#  OutlierDetection.m and OutlierFiltering.m. (See Grané and Veiga, CSDA 2010)
#
#  Input:   alpha   is the significance level,
#           n       is the sample size,
#  Output:  k       is the threshold.
#
#
Threshold=function(alpha,n, shape=7){
    # We generate nMC Monte Carlo samples of size n from the distribution
    # assumed for the model errors. In the following we explicit the case that
    # this assumed distribution is standard normal.

    set.seed(344)
    nMC=20000;
    s=matrix(0,ncol = nMC, nrow = n)

    u=apply(s,2, function(n)(rt(n,shape)))
    # shape=this instruction generates values from a Student's t with 7 degrees of freedom.

    wave <- function(x) {
        wa= dwt(x,filter="haar", n.levels=1)
        w= wa@W[["W1"]]
        sup=max(abs(w))
        return(sup)
    }
    M=apply(u,2,wave)
    k=quantile(M,1-alpha)

    return(k)
}


#####################################################################
#####################################################################
#####################################################################
# Function position=OutlierDetection(u,k)
#
# This function detects the possible outliers (according to a certain
# threshold) in a series of residuals. See Grané and Veiga, CSDA 2010.
#
# Inputs:   u is the series of residuals (in column) to be cleaned,
#           threshold  is the threshold value (see program Threshold.m)
#
# Outputs: position  is a vector containing the positions of the outliers
#                    in the original series.


#12, 125, 456, 667, 799

OutlierDetection <- function(u,k){

    # wavelet decomposition of the residuals
    n=length(u);
    wa= dwt(u,filter="haar", n.levels=1)
    cD=wa@W[["W1"]]
    cA=wa@V[["V1"]]
    l=length(cD); #v=ones(l,1);
    #
    newcD=cD
    JcD=vector()
    coef=vector()
    position=vector()

    # Identification of the outliers in the wavelet coefficients

    j=1
    # Identification of the outliers in the wavelet coefficients
    while(max(abs(newcD)>k)){
        m=max(abs(newcD));
        i=which.max(abs(newcD))
        JcD[j]=i;
        coef[j]=newcD[i];
        newcD[i]=0;
        wa@W[["W1"]]=newcD
        newu = idwt(wa)
        j=j+1;
        wa= dwt(newu,filter="haar", n.levels=1)
        newcD=wa@W[["W1"]]
        cA=wa@V[["V1"]]
    }
    JcD=sort(JcD)
    if(length(JcD) == 0){

        return(list(position=0, JcD=0))}  else { # Identification of the outliers in the original series o

            for(i in 1:length(JcD)){
                u1=na.omit(u[1:(2*abs(JcD[i]-2))])
                u2=na.omit(u[2*(JcD[i]+1):n])
                u3=c(u1,u2)
                mu1=mean(u3)
                if (abs(u[2*JcD[i]]-mu1)>abs(u[(2*JcD[i])-1]-mu1)){
                    position[i]=2*JcD[i];
                }else{
                    position[i]=(2*JcD[i])-1;
                }
            }
            return(list(position=position, JcD=JcD))}
}

#
# Function newseries OutlierFiltering(y,u,k)
#
# This program detects and cleans the outliers in a data series that has been
# previously fitted to a volatility model with completely specified
# distribution for the errors. See Gran? and Veiga, CSDA 2010.
#
#
# Inputs:   y is the data series,
#           u is the series of residuals,
#           k is the threshold value,
# Output:   newseries is the cleaned series.
#
OutlierFiltering <- function(y,u,k){
    n=length(y);
    #----------------------------------------------

    # Decomposition of the series

    wa= dwt(y,filter="haar", n.levels=1)
    cD1=wa@W[["W1"]]
    cA1=wa@V[["V1"]]
    # Identification of the outliers

    wave=OutlierDetection(u,k);

    if (wave$position==0||wave$JcD==0){
        return(list(newseries=y, position=0))
    }else{

        JcD=wave[["JcD"]]
        position=wave[["position"]]
        newcD1=cD1;

        for (i in 1:length(JcD)){
            newcD1[JcD[i]]=0}

        wa@W[["W1"]]=newcD1
        # Reconstruction of the series
        newseries = idwt(wa)
        if (length(position)>1){
            for (l in 1:length(position)){
                if (position[l]<2){
                    newseries[(position[l])+1]=y[(position[l])+1]
                }else if(position[l]<n){
                    newseries[position[l]+1]=y[(position[l]+1)];
                    newseries[position[l]-1]=y[(position[l]-1)];
                }else{
                    newseries[position[l]-1]=y[position[l]-1]
                }
            }
        }
        return(list(newseries=newseries, position=position))
    }
}

############################################################
