
#These functions apply the wavelet-detection when estimating the Value-at-Risk.


VaR_wavelets2<-function(X, res, spec, k){


  new <- OutlierFiltering(X, res, k)

  Y <- new[["newseries"]]


  ctrl = list(n.restarts = 10,tol = 1e-6)

  model1<-try(
    ugarchfit(data = Y, spec = spec, solver = "hybrid",
              solver.control = ctrl),
    silent = TRUE
  )
  if (class(model1) == "try-error"){
    model1 <- ugarchfit(data = Y, spec = spec,
                        solver = "solnp",  #"gosolnp"
                        solver.control = ctrl)
  }
  if (model1@fit$convergence != 0){
    cat("\n\nProblem here!!\n\n")
  }

  Param<-model1@fit[["coef"]]
  spec@model[["fixed.pars"]] <- Param


  H<-ugarchforecast(model1, Y, n.ahead = 1)
  H_1<-as.numeric(H@forecast$sigmaFor)

  ##################################################################################################

  #Wavelets procedue:

  dis = spec@model[["modeldesc"]][["distribution"]]


  if(dis =="norm"){

    Qed = qdist(distribution = dis , mu=0, sigma=1, p=c(0.05, 0.01, 0.005, 0.001))

    f = function(x) (qdist(distribution = dis, mu=0, sigma=1, p = x))

  }else if(dis == "snorm"){

    Qed = qdist(distribution = dis , mu=0, sigma=1, skew = Param["skew"], p=c(0.05, 0.01, 0.005, 0.001))

    f = function(x) (qdist(distribution = dis, mu=0, sigma=1, p = x, skew = Param["skew"]))

  }else if(dis=="std" || dis=="ged"){

    Qed = qdist(distribution = dis , mu=0, sigma=1, shape = Param["shape"], p=c(0.05, 0.01, 0.005, 0.001))

    f = function(x) (qdist(distribution = dis, mu=0, sigma=1, p = x, shape = Param["shape"]))

  }else{

    Qed = qdist(distribution = dis , mu=0, sigma=1, shape = Param["shape"],
                skew = Param["skew"], p=c(0.05, 0.01, 0.005, 0.001))

    f = function(x) (qdist(distribution = dis, mu=0, sigma=1, p = x, skew = Param["skew"], shape = Param["shape"]))

  }


  #Parametric VaR: 5%, 1%, 0.5% and 0.1% level of significance:

  # Considering the shape parameter estimated previously.

  VaR <- as.numeric(H_1*(Qed))

  #######################################################################
  #Expected shortfall (ES)
  #######################################################################


  ES5a =  H_1*integrate(f, 0, 0.05)$value/0.05

  ES1a =  H_1*integrate(f, 0, 0.01)$value/0.01

  ES5t =  H_1*integrate(f, 0, 0.005)$value/0.005

  ES1t =  H_1*integrate(f, 0, 0.001)$value/0.001



  ans<-list(

    VaR5a=VaR[1],
    VaR1a=VaR[2],
    VaR5t=VaR[3],
    VaR1t=VaR[4],

    #VaR5F=as.numeric(FHS1[1]),
    #VaR1F=as.numeric(FHS1[2]),
    #VaR5tF=as.numeric(FHS1[3]),
    #VaR1tF=as.numeric(FHS1[4]),

    ES5a = as.numeric(ES5a),
    ES1a = as.numeric(ES1a),
    ES5t = as.numeric(ES5t),
    ES1t = as.numeric(ES1t),

    Vol_w = H_1


  )

  return(ans)

}

######################################################################################3
#######################################################################################
#######################################################################################


VaR_wavelets3 <- function(X, res, model, k){


  new <- OutlierFiltering(X, res, k)

  Y <- new[["newseries"]]


  if (model == "NGQML"){

    model1 <- NGQML_f(Y)

   }else if (model == "Rank"){

    model1 <- Rank_f(Y)

   }else if (model == "RankB"){

     model1 <- Rank_B(Y)

   }else if (model == "Prem") {

     model1 <- Prem_f1(Y)

  }

  specw = ugarchspec(variance.model = list(model = "sGARCH"),
                    mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
                    distribution.model = dis)

  specw@model[["fixed.pars"]]<- model1


  H<-ugarchforecast(specw, Y, n.ahead = 1)
  H_1<-as.numeric(H@forecast$sigmaFor)

  ##################################################################################################
  Qed = qdist(distribution = dis , mu=0, sigma=1, shape = 7,
              p=c(0.05, 0.01, 0.005, 0.001))

  f = function(x) (qdist(distribution = dis, mu=0, sigma=1, p = x, shape = 7))

  ES5a = H_1*integrate(f, 0, 0.05)$value/0.05

  ES1a = H_1*integrate(f, 0, 0.01)$value/0.01

  ES5t = H_1*integrate(f, 0, 0.005)$value/0.005

  ES1t = H_1*integrate(f, 0, 0.001)$value/0.001

  VaR <- as.numeric(H_1*(Qed))

  #################################################################################

  ## Wavelets



  ans<-list(

    VaR5w = VaR[1],
    VaR1w = VaR[2],
    VaR5tw = VaR[3],
    VaR1tw = VaR[4],


    ES5w = as.numeric(ES5a),
    ES1w = as.numeric(ES1a),
    ES5tw = as.numeric(ES5t),
    ES1tw = as.numeric(ES1t),


    Volw = H_1


  )
}
