
#main author: Kévin Allan Sales Rodrigues

#### algorithm to make inference in LAD models

#' Fitting LAD Models
#'
#' @param y A vector with response variables.
#' @param X A matrix or vector with explanatory variables.
#' @param intercept 1 for a model with intercept and 0 for a model without intercept.
#' @param alpha  significance level to hypotheses tests.
#' @param print 1 to print response variables, fitted values and residuals, 0 to don't print.
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
#' ladreg(stack.loss, stack.x, intercept =1, alpha=0.05, print=1)


#Variáveis a serem definidas pelo usuário:

#intercept: recebendo o valor 0, ajusta o modelo sem intercepto e recebendo o valor 1, ajusta o modelo com intercepto.

#print: recebendo o valor 1, imprime os valores observados, ajustados e os resíduos do modelo e recebendo 0, não imprime tais valores.



#Verifica se o usuario pediu modelo com ou sem intercepto

##colocando valores "default"
ladreg = function(y,X, intercept =1, alpha=0.05, print=1){

  if(intercept == 0){
    fit = ladfit(X, y, intercept = F)
    Xmat =X
    SEAm = sum(abs(y))
    fator=0
  }

  if(intercept == 1){
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
  IC[,1] = coef1 -qnorm(0.975)*tao*aux
  IC[,2] = coef1 +qnorm(0.975)*tao*aux
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
    if(intercept == 0) (fit1 =ladfit(X1,y1, intercept = F))
    if(intercept == 1) (fit1 =ladfit(X1,y1, intercept = T))
    coefaux = cbind(coef(fit1))
    yaju = Xmat[i,]%*%coefaux
    parc = abs(y[i] -yaju)
    SEAP = SEAP +parc
  }

  #Observacoes Definidoras
  print(def)

  #Verifica se o usuario pediu a impressao dos valores observados, ajustados e residuos
  if(print == 1) {
    OAR = cbind(y, aju, r)
    #Valores Observados, Ajustados e Residuos
    print(OAR)
  }

  #Estimativa de tao
  tao = round(tao, digits=4)

  #Estimativas dos Parametros do Modelo
  coef1 = round(coef1, digits=4)
  print(coef1)

  #Estimativas dos Erros Padrao dos Estimadores do Modelo
  DP = round(DP, digits=4)
  print(DP)

  #Intervalos de Confianca (95%)
  IC = round(IC, digits =4)
  print(IC)

  #Estatistica do Teste Individual dos Parametros e Nivel Descritivo
  TH = round(TH, digits =4)
  print(TH)

  #Estatistica do Teste do Efeito de Regressao e Nivel Descritivo
  EReg = round(EReg, digits =4)
  print(EReg)

  ###Estatistica para Selecao de Variaveis

  #Soma dos Erros Absolutos (SEA)
  SEA = round(SEA, digits=4)

  #Erro Absoluto Medio (EAM)
  EAM = round(EAM, digits =4)

  #Coeficiente de Determinacao Multipla - R2
  R2 = round(R2, digits=4)

  #Soma do Erros absolutos Preditivos (SEAP)
  SEAP =round(SEAP, digits=4)

  FIT = list("OAR" = OAR, "tao" = tao, "coef" = coef1,
             "DP" = DP, "IC" = IC,  "TH" = TH, "EReg" = EReg,
             "SEA" = SEA, "EAM" = EAM, "R2" = R2, "SEAP" = SEAP)

  return(FIT)

}
