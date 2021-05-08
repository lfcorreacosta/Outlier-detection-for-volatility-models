library(xts)
library(fBasics)


tables_paper <- function(matr1, matr2, stati){

FZa = matr1[[3]][[1]]
magi = matr1[[1]]

upa0 = cbind(magi, FZa)
upa = t(upa0)

i = 1; j =2; ki = 3; l = 4; m = 5; n = 6
stati = stati

models <- c("sGARCH","iGARCH", "eGARCH", "gjrGARCH", "apARCH", "csGARCH")
modelsW <- paste(models, "wav", sep="-")
distributions <- c("norm", "std", "ged", "snorm", "sstd", "sged")

LRucA<-cbind( rbind(upa[stati,i], upa[stati,(i+6)], upa[stati,(i+6*2)], upa[stati,(i+6*3)], upa[stati,(i+6*4)], upa[stati,(i+6*5)]),
rbind(upa[stati,j], upa[stati,(j+6)], upa[stati,(j+6*2)], upa[stati,(j+6*3)], upa[stati,(j+6*4)], upa[stati,(j+6*5)]),
LRucG = rbind(upa[stati,ki], upa[stati,(ki+6)], upa[stati,(ki+6*2)], upa[stati,(ki+6*3)], upa[stati,(ki+6*4)], upa[stati,(ki+6*5)]),
rbind(upa[stati,l], upa[stati,(l+6)], upa[stati,(l+6*2)], upa[stati,(l+6*3)], upa[stati,(l+6*4)], upa[stati,(l+6*5)]),
rbind(upa[stati,m], upa[stati,(m+6)], upa[stati,(m+6*2)], upa[stati,(m+6*3)], upa[stati,(m+6*4)], upa[stati,(m+6*5)]),
rbind(upa[stati,n], upa[stati,(n+6)], upa[stati,(n+6*2)], upa[stati,(n+6*3)], upa[stati,(n+6*4)], upa[stati,(n+6*5)]))

colnames(LRucA) = distributions
rownames(LRucA) = models
 

FZw = matr2[[3]][[1]]
magw = matr2[[1]]

upw0 = cbind(magw, FZw)
upw = t(upw0)


LRucW<-cbind( rbind(upw[stati,i], upw[stati,(i+6)], upw[stati,(i+6*2)], upw[stati,(i+6*3)], upw[stati,(i+6*4)], upw[stati,(i+6*5)]),
              rbind(upw[stati,j], upw[stati,(j+6)], upw[stati,(j+6*2)], upw[stati,(j+6*3)], upw[stati,(j+6*4)], upw[stati,(j+6*5)]),
              LRucG = rbind(upw[stati,ki], upw[stati,(ki+6)], upw[stati,(ki+6*2)], upw[stati,(ki+6*3)], upw[stati,(ki+6*4)], upw[stati,(ki+6*5)]),
              rbind(upw[stati,l], upw[stati,(l+6)], upw[stati,(l+6*2)], upw[stati,(l+6*3)], upw[stati,(l+6*4)], upw[stati,(l+6*5)]),
              rbind(upw[stati,m], upw[stati,(m+6)], upw[stati,(m+6*2)], upw[stati,(m+6*3)], upw[stati,(m+6*4)], upw[stati,(m+6*5)]),
              rbind(upw[stati,n], upw[stati,(n+6)], upw[stati,(n+6*2)], upw[stati,(n+6*3)], upw[stati,(n+6*4)], upw[stati,(n+6*5)]))

colnames(LRucW) = distributions
rownames(LRucW) = modelsW


Teste3 <- cbind(as.numeric(matr1[[1]][,4]), as.numeric(matr1[[3]][[1]]), as.numeric(matr1[[1]][,c(5,6)]),
                as.numeric(matr2[[1]][,4]), as.numeric(matr2[[3]][[1]]), as.numeric(matr2[[1]][,c(5,6)]))

Teste3 <- cbind(names(VaR1a),matr1[[1]][,4], matr1[[3]][[1]], matr1[[1]][,c(5,6)],
                names(VaR1w),matr2[[1]][,4], matr2[[3]][[1]], matr2[[1]][,c(5,6)])

colnames(Teste3) <-rep(c("models", "A/E","FZ", "ADMin", "ADMax"),2)
rownames(Teste3)=NULL


ans = rbind(LRucA, LRucW)

return(list(test=ans, lossF = Teste3))

}

