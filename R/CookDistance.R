#main author: Kévin Allan Sales Rodrigues

#### uma medida de influencia

#' Calculate Cook Distance
#'
#' @param X A matrix or vector with explanatory variables.
#' @param y A vector with response variables.
#' @return
#' \item{Cook Distance}{          A vector with Cook Distance for each observation.}.
#'
#' @references Sun, R.-B. and Wei, B.-C. (2004) On influence assessment for lad regression.
#' \emph{Statistics & Probability Letters}, \strong{67}, 97-110. \doi{10.1016/j.spl.2003.08.018}.
#'
#'
#' @examples
#' ### Using stackloss data
#'
#' CookDistance(stack.loss, stack.x)


CookDistance = function(y,X){

  n = length(y)

  ###Calculando taoi que e o tao sem a i-esima observacao
  ###e calculando os betai que sao betai[,i]

  beta = coef(ladfit(as.matrix(X), y))

  betai = coef(ladfit(as.matrix(X)[-1,], y[-1]))
  taoi = sum(abs(y-cbind(1,X)%*%as.matrix(betai)[,1]))

  for (i in 2:n){
    betai = cbind(betai, coef(ladfit(as.matrix(X)[-i,], y[-i])))
    taoi = c(taoi, sum(abs(y-cbind(1,X)%*%as.matrix(betai)[,i])))
  }
  taoi = taoi/(n-1)

  ###Calculando LD(betai, taoi)

  ##por enquanto
  #tao = 1.5088

  tao = sum(abs(y-cbind(1,X)%*%as.matrix(coef(ladfit(X, y)))[,1]))/n


  CD = tao^(-2)* t(beta-betai[,1])%*%(t(cbind(1,X)) %*%cbind(1,X))%*% (beta-betai[,1])

  for (i in 2:n){
    CD = c(CD, tao^(-2)* t(beta-betai[,i])%*%(t(cbind(1,X)) %*%cbind(1,X))%*% (beta-betai[,i]))
  }

  plot(CD, main="Cook Distance", xlab="Indices",ylab="Cook Distance")
  return(CD)
  #return(beta%*%t(rep(1,n))-betai)
  #return((t(cbind(1,X)) %*%cbind(1,X)))
}

