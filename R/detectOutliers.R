#main author: Kévin Allan Sales Rodrigues

#### algorithm to detect outliers

#' Algorithm To Detect Outliers
#'
#' @param X A matrix or vector with explanatory variables.
#' @param y A vector with response variables.
#' @param intercepto 1 for a model with intercept and 0 for a model without intercept.
#' @return
#' \item{outliers     }{          A vector with outliers indices.}.
#'
#' @references Dodge, Y. (1997) Lad regression for detecting outliers in response and explanatory variables.
#' \emph{Journal of Multivariate Analysis}, \strong{61}, 144–158. \doi{10.1006/jmva.1997.1666}
#'
#'
#'
#' @examples
#' ### Using stackloss data
#'
#' detectOutliers(stack.loss, stack.x)

##colocando valores "default"
detectOutliers= function(y,X, intercepto =1){

residuos =  ladfit(X,y)$residuals
#print(ladfit(X,y)$coef)

##número de observações
n = length(y)
##numero de parametros
p = ncol(as.matrix(X)) + intercepto

###ordenando resíduos absolutos em ordem decrescente
orde =sort(abs(residuos),decreasing = TRUE)

###excluindo resíduos nulos
result =NULL
for(i in 1:(n-p)){
result[i]=orde[i]
}

##MRADZ é a mediana dos resíduos absolutos diferentes de zero
##estimativa do lambda
lambda = 1.4826*stats::median(result)

#print(lambda)

##calculando resíduos "padronizados"
residuos2 = residuos/lambda


##checando se os resíduos "padronizados" > 2,5
residuos2 = residuos/lambda
#print(residuos2)

##criando conjunto com as observações excluidas da análise
excluir =c()

for(i in 1:length(residuos2)){
  if(abs(residuos2[i]) >= 2.5){
    excluir = c(excluir,i)
  }
}
  ###etapa 1 encerrada (excluindo outliers em Y)

###etapa 2 (excluindo outliers em X)
if(!is.null(excluir)){
  n = length(y[-excluir])
}else{
  n = length(y)
}

##criando conjunto com as observações excluidas da análise
excluir2 =c()

for(i in 1:ncol(as.matrix(X))){

  if(!is.null(excluir)){
    residuos =  ladfit(cbind(y[-excluir], as.matrix(X)[-excluir,-i]),as.matrix(X)[-excluir,i])$residuals
  }else{
    residuos =  ladfit(cbind(y, as.matrix(X)[,-i]),as.matrix(X)[,i])$residuals
  }

#  if(!is.null(excluir)){
#    print(ladfit(cbind(y[-excluir], as.matrix(X)[-excluir,-i]),as.matrix(X)[-excluir,i])$coef)
#  }else{
#    print(ladfit(cbind(y, as.matrix(X)[,-i]),as.matrix(X)[,i])$coef)
#  }


  ###ordenando resíduos absolutos em ordem decrescente
  orde =sort(abs(residuos),decreasing = TRUE)

  ###excluindo resíduos nulos
  result =NULL
  for(j in 1:(n-p)){
    result[j]=orde[j]
  }

  ##MRADZ é a mediana dos resíduos absolutos diferentes de zero
  ##estimativa do lambda
  lambda = 1.4826*stats::median(result)

#print(lambda)

  ##calculando resíduos "padronizados"
  residuos2 = residuos/lambda


  ##checando se os resíduos "padronizados" > 2,5
  residuos2 = residuos/lambda
  #print(residuos2)

  #####criando conjunto com as observações excluidas da análise (outliers em X)

  for(k in 1:length(residuos2)){
    if(abs(residuos2[k]) >= 2.5){
      excluir2 = c(excluir2,k)
    }
  }
  #print(excluir2)

} ### fim da etapa 2 (excluindo outliers em X)


##checando quais observações restaram no modelo
indices = 1:length(y)

##testando se o conjunto não é nulo
if(!is.null(excluir)){
  indices =indices[-excluir]
}

##testando se o conjunto não é nulo
if(!is.null(excluir2)){
  excluir2 =unique(excluir2)
  indices =indices[-unique(excluir2)]
}

###determinando outilers
outliers=setdiff(1:length(y),indices)


###esse vetor contém as observações que restaram
return(outliers)
}


