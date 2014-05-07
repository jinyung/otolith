# aligne configuration
# 
# @description aligne the landmark configuration to the PC axis
# @details modified from J Claude's Book Function, which in the book, 
#  the \code{aligne} function takes array, here, the \code{aligne2} function
# takes matrix
# @references Claude J. (2008). Morphometrics with R. Springer

aligne2<-function(A){
  B<-A
  Ms<-scale(A, scale=F)
  sv<-eigen(var(Ms))
  M<-Ms%*%sv$vectors
  B<-M
  B
} 