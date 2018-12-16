#main author: Kévin Allan Sales Rodrigues

#### medida de influencia

#' Calculate Likelihood Displacement
#'
#' @param X A matrix or vector with explanatory variables.
#' @param y A vector with response variables.
#' @return
#' \item{Likelihood Displacement}{          A vector with Likelihood Displacement for each observation.}.
#'
#' @references Elian, S. N., André, C. D. S. and Narula, S. C. (2000) Influence Measure for the
#' \ifelse{html}{\out{L<sub>1</sub>}}{\eqn{L1}} regression.
#' \emph{Communications in Statistics - Theory and Methods}, \strong{29}(4), 837-849. \doi{10.1080/03610920008832518}.
#'
#'
#' @examples
#' ### Using stackloss data
#'
#' likelihoodD(stack.loss, stack.x)


likelihoodD = function(y,X){

  n = length(y)

  ###Calculando taoi que e o tao sem a i-esima observacao
  ###e calculando os betai que sao betai[,i]


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

  LD = 2*(n*log(taoi[1]/tao) + abs(y[1] -t(as.matrix(cbind(1,X)[1,]))%*%as.matrix(betai)[,1])/taoi[1] -1)

  for (i in 2:n){
    LD = c(LD, 2*(n*log(taoi[i]/tao) + abs(y[i] -t(as.matrix(cbind(1,X)[i,]))%*%as.matrix(betai)[,i])/taoi[i] -1))
  }

  plot(LD, main="Likelihood Displacement", xlab="Indices",ylab="Likelihood Displacement")
  return(LD)
}


