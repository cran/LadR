#main author: Kévin Allan Sales Rodrigues

#### medida de influencia condicional

#' Calculate Conditional Likelihood Displacement
#'
#' @param X A matrix or vector with explanatory variables.
#' @param y A vector with response variables.
#' @return
#' \item{Conditional Likelihood Displacement}{          A vector with Conditional Likelihood Displacement for each observation.}.
#'
#' @references Elian, S. N., André, C. D. S. and Narula, S. C. (2000) Influence Measure for the
#' \ifelse{html}{\out{L<sub>1</sub>}}{\eqn{L1}} regression.
#' \emph{Communications in Statistics - Theory and Methods}, \strong{29}(4), 837-849. \doi{10.1080/03610920008832518}.
#'
#'
#' @examples
#' ### Using stackloss data
#'
#' likelihoodDC(stack.loss, stack.x)


likelihoodDC = function(y,X){

  n = length(y)

  betai = coef(ladfit(as.matrix(X)[-1,], y[-1]))
  taoi = sum(abs(y-cbind(1,X)%*%as.matrix(betai)[,1]))

  for (i in 2:n){
    betai = cbind(betai, coef(ladfit(as.matrix(X)[-i,], y[-i])))
    taoi = c(taoi, sum(abs(y-cbind(1,X)%*%as.matrix(betai)[,i])))
  }
  taoi = taoi/(n-1)

  ###Calculando LD(betai| taoi)
  beta = coef(ladfit(X, y))

  ##por enquanto
  #tao = 1.5088

  tao = sum(abs(y-cbind(1,X)%*%as.matrix(coef(ladfit(X, y)))[,1]))/n

  LD_condicional =  2*n*log(sum(abs(y -as.matrix(cbind(1,X))%*%as.matrix(betai)[,1]))/sum(abs(y -as.matrix(cbind(1,X))%*%as.matrix(beta))))

  for (i in 2:n){
    LD_condicional = c(LD_condicional, 2*n*log(sum(abs(y -as.matrix(cbind(1,X))%*%as.matrix(betai)[,i]))/sum(abs(y -as.matrix(cbind(1,X))%*%as.matrix(beta)))))
  }

  plot(LD_condicional, main="Likelihood Displacement Condicional", xlab="Indices",ylab="Likelihood Displacement")
  return(LD_condicional)

}

