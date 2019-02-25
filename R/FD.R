#main author: Kévin Allan Sales Rodrigues

#### uma medida de influencia

#' F Distance
#'
#' @param X A matrix or vector with explanatory variables.
#' @param y A vector with response variables.
#' @return
#' \item{F Distance}{          A vector with F Distance for each observation.}.
#'
#' @references Sun, R.-B. and Wei, B.-C. (2004) On influence assessment for lad regression.
#' \emph{Statistics & Probability Letters}, \strong{67}, 97-110. \doi{10.1016/j.spl.2003.08.018}.
#'
#'
#' @examples
#' ### Using stackloss data
#'
#' CookDistance(stack.loss, stack.x)


FD = function(y,X){

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


  FF = tao^(-2)*sqrt(sum((cbind(1,X)%*%beta -cbind(1,X)%*%betai[,1])^(2)))

  for (i in 2:n){
    FF = c(FF, tao^(-2)*sqrt(sum((cbind(1,X)%*%beta -cbind(1,X)%*%betai[,i])^(2))))
  }


  plot(FF, main="F Distance", xlab="Indices",ylab="F Distance")
  return(FF)
}

