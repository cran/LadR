
#main author: Kévin Allan Sales Rodrigues

#### algorithm to make inference in LAD models

#' Fitting LAD Models
#'
#' @param y A vector with response variables.
#' @param X A matrix or vector with explanatory variables.
#' @param intercepto 1 for a model with intercept and 0 for a model without intercept.
#' @param alfa  significance level to hypotheses tests.
#' @param imprime 1 to print response variables, fitted values and residuals, 0 to don't print.
#' @return
#' \item{outliers     }{          A vector with outliers indices.}.
#'
#' @references Dielman, T. E. (2005) Least absolute value regression: recent contributions.
#' \emph{Journal of Statistical Computational and Simulation}, \strong{75}(4), 263–286. \doi{10.1080/0094965042000223680}
#'
#'
#'
#' @examples
#' ### Using stackloss data
#'
#' ladreg(stack.loss, stack.x, intercepto =1, alfa=0.05, imprime=1)


#Variáveis a serem definidas pelo usuário:

#intercepto: recebendo o valor 0, ajusta o modelo sem intercepto e recebendo o valor 1, ajusta o modelo com intercepto.

#imprime: recebendo o valor 1, imprime os valores observados, ajustados e os resíduos do modelo e recebendo 0, não imprime tais valores.

#laplace: recebendo o valor 1, imprime as estatísticas de seleção de variáveis que podem ser usadas no caso da distribuição dos dados
#ser de Laplace e recebendo 0, não imprime tais estatísticas.


#Verifica se o usuario pediu modelo com ou sem intercepto

##colocando valores "default"
ladreg = function(y,X, intercepto =1, alfa=0.05, imprime=1){

  laplace =0

  if(intercepto == 0){
    fit = ladfit(X, y, intercept = F)
    Xmat =X
    SEAm = sum(abs(y))
    fator=0
  }

  if(intercepto == 1){
    fit = ladfit(X, y, intercept = T)
    Xmat =cbind(1,X)
    SEAm = sum(abs(y-median(y)))
    fator=1
  }

  #Calculo da Estimativa de tao
  r= cbind(residuals(fit))
  coef1 = cbind(coef(fit))
  aju = y-r
  n = nrow(X) - ncol(X) -fator
  m = (n+1)/2 -sqrt(n)
  res = cbind(sort(r))
  res1 = matrix(0, nrow(res) -nrow(coef1),1)
  def = matrix(0, nrow(coef1),1)
  j=1
  for(i in 1:nrow(r)) if(r[i] == 0) ((def[j] = i) & (j = j+1))
  j=1
  for(i in 1:nrow(res)) if(res[i] !=0) ((res1[j] = res[i]) & (j = j+1))
  pos1 = n-m +1
  pos2 = m
  e1 = round(pos1)
  e2 = round(pos2)
  if (e1 > n) (e1 = n)
  if(e2 == 0) (e2 = 1)
  tao = sqrt(n)*(res[e1] -res1[e2])/4

  #Calculos dos Intervalos de Confianca, Estatisticas de Testes e Niveis descritivos

  SEA = sum(abs(res))
  XX = solve(t(Xmat)%*%Xmat)
  DX = diag(XX)
  IC = matrix(0, nrow(coef1),2)
  TH = matrix(0, nrow(coef1),2)
  #DP = matrix(0, nrow(coef1),1)
  aux = cbind(sqrt(DX))
  DP = tao*aux
  IC[,1] = coef1 -1.96*tao*aux
  IC[,2] = coef1 +1.96*tao*aux
  E = matrix(0, nrow(coef1), 1)
  for(i in 1:ncol(Xmat)) E[i] = coef1[i]/(tao*aux[i])
  E = abs(E)
  TH[,1] = E
  TH[,2] = 2*(1-pnorm(E,0,1))
  EReg = matrix(0,1,2)
  reg = (SEAm -SEA)/(tao/2)
  EReg[1,1] = reg
  EReg[1,2] = 1-stats::pchisq(reg, ncol(X))

  #Calculos das Estatisticas para Selecao de Variaveis

  EAM = SEA/length(y)
  RSEA = SEAm -SEA
  R2 = RSEA/(RSEA +(nrow(X) -ncol(X) -fator)*tao/2)
  R3 = RSEA/(RSEA +(nrow(X) -ncol(X) -fator)*EAM/2)
  Raju = 1 -(1-R3)*(nrow(y)/(nrow(X) -ncol(X) -fator))
  SEAP = 0
  for(i in 1:nrow(X)) {
    X1 = matrix(0, nrow(X) -1, ncol(X))
    y1 = matrix(0, nrow(X) -1, 1)
    k=1
    for(j in 1:nrow(X)) if (j != i) (X1[k,] = X[j,]) & (y1[k] = y[j]) & (k=k+1)
    if(intercepto == 0) (fit1 =ladfit(X1,y1, intercept = F))
    if(intercepto == 1) (fit1 =ladfit(X1,y1, intercept = T))
    coefaux = cbind(coef(fit1))
    yaju = Xmat[i,]%*%coefaux
    parc = abs(y[i] -yaju)
    SEAP = SEAP +parc
  }

  #Impressao dos Resultados
  cat(" ", fill=T)
  cat("Observacoes Definidoras:", fill=T)
  print(def)
  cat(" ", fill=T)

  #Verifica se o usuario pediu a impressao dos valores observados, ajustados e residuos
  if(imprime == 1) {
    OAR = cbind(y, aju, r)
    cat("Valores Observados, Ajustados e Residuos:", fill =T)
    cat(" ")
    print(OAR)
    cat(" ", fill=T)
  }
  tao = round(tao, digits=4)
  cat("Estimativa de tao:", tao, fill=T)
  cat(" ", fill=T)
  cat("Estimativas dos Parametros do Modelo:", fill=T)
  coef1 = round(coef1, digits=4)
  print(coef1)
  cat(" ", fill=T)
  cat("Estimativas dos Erros Padrao dos Estimadores do Modelo:", fill=T)
  DP = round(DP, digits=4)
  print(DP)
  cat(" ", fill=T)
  cat("Intervalos de Confianca (95%):", fill=T)
  IC = round(IC, digits =4)
  print(IC)
  cat(" ", fill=T)
  cat("Estatistica do Teste Individual dos Parametros e Nivel Descritivo:", fill=T)
  TH = round(TH, digits =4)
  print(TH)
  cat(" ", fill=T)
  cat("Estatistica do Teste do Efeito de Regressao e Nivel Descritivo:", fill=T)
  EReg = round(EReg, digits =4)
  print(EReg)
  cat(" ", fill=T)
  cat(" ", fill=T)
  cat("Estatistica para Selecao de Variaveis:", fill=T)
  cat(" ", fill=T)
  SEA = round(SEA, digits=4)
  cat("Soma dos Erros Absolutos (SEA):", SEA, fill=T)
  cat(" ", fill=T)
  EAM = round(EAM, digits =4)
  cat("Erro Absoluto Medio (EAM):", EAM, fill=T)
  cat(" ", fill=T)
  R2 = round(R2, digits=4)
  cat("Coeficiente de Determinacao Multipla - R2:", R2, fill=T)
  cat(" ", fill=T)

  #Verifica se o usuario pediu a impressao das Estatisticas usadas para a Distribuicao Laplace

  if(laplace == 1){
    R3 = round(R3, digits=4)
    cat("Coeficiente de Determinacao Multipla - R3:", R3, fill=T)
    cat(" ", fill=T)
    Raju = round(Raju, digits=4)
    cat("Coeficiente de Determinacao Multipla Ajustado:", Raju, fill=T)
    cat(" ", fill=T)
  }

  SEAP =round(SEAP, digits=4)
  cat("Soma do Erros absolutos Preditivos (SEAP):", SEAP, fill=T)
  cat(" ", fill=T)
  cat("**********************************************************", fill=T)
  cat(" ", fill=T)

}
