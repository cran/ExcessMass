lambdaweight <-
function(xdata,lambda){
  nlambda <- length(xdata)*lambda
  nHLambda <- (1:length(xdata))-xdata*nlambda
  nHLambda
}
