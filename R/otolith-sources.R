# ====== FUNCTIONS SOURCED FROM OTHER PLACES ====== #
# 1. JClaude Book functions                  
# 2. kfcv.sizes() and kfcv.testing()         
# 3. aligne2()                                

#-------------------------------------------------------------
# JClaude Book functions
# @author Julien Claude
# @references Claude J. (2008). Morphometrics with R. Springer
# @description functions developped in the book.
# @note tps function is modified
#-------------------------------------------------------------

######functions developped in Claude J. (2008). Morphometrics with R. Springer. 
######Refer to the book if you use them
####f2.8
eigenrotation<-function(N, M)
{M<-as.matrix(M); N<-as.matrix(N)
 sN<-eigen(var(N))$vectors
 k<-dim(N)[2]
 p<-dim(M)[1]
 Nn<-N%*%sN
 Mn<-M%*%sN
 uNn<-apply(Nn, 2, mean)
 Mnf<-Mn-rep(1,p)%*%t(uNn)
 Mnf}
####f4.1
centcoord<-function(M){apply(M,2,mean)}
####f4.2
centsiz<-function(M)
{p<-dim(M)[1]
 size<-sqrt(sum(apply(M, 2,var))*(p-1))
 list("centroid_size" = size,"scaled" = M/size)}
####f4.10
transl<-function(M)
{M - matrix(centcoord(M), nrow(M), ncol(M), byrow=T)}
####f4.11
trans1<-function(M){scale(M,scale=F)}
####f4.12
helmert<-function(p)
{H<-matrix(0, p, p)
 diag(H)<--(0:(p-1)) * (-((0:(p-1))*((0:(p-1))+1))^(-0.5))
 for (i in 2:p){H[i,1:(i-1)]<- -((i-1)*(i))^(-0.5)}
 H[1,]<-1/sqrt(p)
 H}
####f4.13
helmertm<-function(M)
{helmert(nrow(M))[-1,]%*%M}
####f4.14
kendall2d<-function(M)
{Mh<-helmertm(M)
 mhc<-complex(nrow(Mh),Mh[,1], Mh[,2])
 (mhc/mhc[1])[-1]}
####f4.15
fPsup<-function(M1, M2)
{k<-ncol(M1)
 Z1<-trans1(centsiz(M1)[[2]])
 Z2<-trans1(centsiz(M2)[[2]])
 sv<-svd(t(Z2)%*%Z1)
 U<-sv$v; V<-sv$u; Delt<-sv$d
 sig<-sign(det(t(Z2)%*%Z1))
 Delt[k]<-sig*abs(Delt[k]) ; V[,k]<-sig * V[,k]
 Gam<-U%*%t(V)
 beta<-sum(Delt)
 list(Mp1=beta*Z1%*%Gam,Mp2=Z2,rotation=Gam,scale=beta,DF=sqrt(1-beta^2))}
####f4.16
ild2<-function(M1, M2){sqrt(apply((M1-M2)^2, 1, sum))}
####f4.17
pPsup<-function(M1,M2)
{k<-ncol(M1)
 Z1<-trans1(centsiz(M1)[[2]])
 Z2<-trans1(centsiz(M2)[[2]])
 sv<-svd(t(Z2)%*%Z1)
 U<-sv$v; V<-sv$u; Delt<-sv$d
 sig<-sign(det(t(Z1)%*%Z2))
 Delt[k]<-sig*abs(Delt[k]) ; V[,k]<-sig * V[,k]
 Gam<-U%*%t(V)
 beta<-sum(Delt)
 list(Mp1=Z1%*%Gam,Mp2=Z2, rotation=Gam,DP=sqrt(sum(ild2(Z1%*%Gam, Z2)^2)),rho=acos(beta))}
####f4.18
mshape<-function(A){apply(A, c(1,2), mean)}
####f4.19
fgpa<-function(A){
  p<-dim(A)[1]; k<-dim(A)[2]; n<-dim(A)[3]
  temp2<-temp1<-array(NA, dim=c(p,k,n))
  Siz<-numeric(n)
  for(i in 1:n)
  {Acs<-centsiz(A[,,i])
   Siz[i]<-Acs[[1]]
   temp1[,,i]<-trans1(Acs[[2]])}
  Qm1<-dist(t(matrix(temp1,k*p,n)))
  Q<-sum(Qm1); iter<-0
  while (abs(Q)>0.00001)
  {for (i in 1:n){
    M<-mshape(temp1[,,-i])
    temp2[,,i]<-fPsup(temp1[,,i],M)[[1]]}
   Qm2<-dist(t(matrix(temp2,k*p,n)))
   Q<-sum(Qm1)-sum(Qm2)
   Qm1<-Qm2
   iter=iter+1
   temp1<-temp2}
  list(rotated=temp2,iterationnumber=iter,Q=Q,interproc.dist=Qm2,mshape=centsiz(mshape(temp2))[[2]],cent.size=Siz)}
####f4.20
fgpa2<-function(A)
{p<-dim(A)[1]; k<-dim(A)[2]; n<-dim(A)[3]
 temp2<-temp1<-array(NA, dim=c(p,k,n))
 Siz<-numeric(n)
 for (i in 1:n)
 {Acs<-centsiz(A[,,i])
  Siz[i]<-Acs[[1]]
  temp1[,,i]<-trans1(Acs[[2]])}
 iter<-0; sf<-NA
 M<-temp1[,,1]
 for (i in 1:n)
 {temp1[,,i]<-fPsup(temp1[,,i],M)[[1]]}
 M<-mshape(temp1)
 Qm1<-dist(t(matrix(temp1,k*p,n)))
 Q<-sum(Qm1); iter<-0
 sc<-rep(1,n)
 while (abs(Q)>0.00001){
   for (i in 1:n){
     Z1<-temp1[,,i]
     sv<-svd(t(M)%*%Z1)
     U<-sv$v; V<-sv$u; Delt<-sv$d
     sig<-sign(det(t(Z1)%*%M))
     Delt[k]<-sig*abs(Delt[k])
     V[,k]<-sig*V[,k]
     phi<-U%*%t(V)
     beta<-sum(Delt)
     temp1[,,i]<-X<-sc[i]*Z1%*%phi}
   M<-mshape(temp1)
   for (i in 1:n)
   {sf[i]<-sqrt(sum(diag(temp1[,,i]%*%t(M)))
                /(sum(diag(M%*%t(M)))*sum(diag(temp1[,,i]
                                               %*%t(temp1[,,i])))))
    temp2[,,i]<-sf[i]*temp1[,,i]}
   M<-mshape(temp2)
   sc<-sf*sc
   Qm2<-dist(t(matrix(temp2,k*p,n)))
   Q<-sum(Qm1)-sum(Qm2)
   Qm1<-Qm2
   iter=iter+1
   temp1<-temp2}
 list(rotated=temp2,iterationnumber=iter,Q=Q,intereuclidean.dist=Qm2, mshape=centsiz(mshape(temp2))[[2]], cent.size=Siz)}
####f4.21
pgpa<-function(A)
{p<-dim(A)[1];k<-dim(A)[2];n<-dim(A)[3]
 temp2<-temp1<-array(NA, dim=c(p,k,n)); Siz<-numeric(n)
 for (i in 1:n)
 {Acs<-centsiz(A[,,i])
  Siz[i]<-Acs[[1]]
  temp1[,,i]<-trans1(Acs[[2]])}
 Qm1<-dist(t(matrix(temp1,k*p,n)))
 Q<-sum(Qm1); iter<-0
 while (abs(Q)>0.00001)
 {for(i in 1:n){
   M<-mshape(temp1[,,-i])
   temp2[,,i]<-pPsup(temp1[,,i],M)[[1]]}
  Qm2<-dist(t(matrix(temp2,k*p,n)))
  Q<-sum(Qm1)-sum(Qm2)
  Qm1<-Qm2
  iter=iter+1
  temp1<-temp2}
 list("rotated"=temp2,"it.number"=iter,"Q"=Q,"intereucl.dist"=Qm2,"mshape"=centsiz(mshape(temp2))[[2]],"cent.size"=Siz)}
####f4.22#note the function aligne is the book is bugged. This one is working.
aligne<-function(A)
{B<-A
 n<-dim(A)[3]; k<-dim(A)[2]
 for (i in 1:n)
 {Ms<-scale(A[,,i], scale=F)
  sv<-eigen(var(Ms))
  M<-Ms%*%sv$vectors
  B[,,i]<-M}
 B}
####f4.24
orp<-function(A)
{p<-dim(A)[1];k<-dim(A)[2];n<-dim(A)[3]
 Y1<-as.vector(centsiz(mshape(A))[[2]])
 oo<-as.matrix(rep(1,n))%*%Y1
 I<-diag(1,k*p)
 mat<-matrix(NA, n, k*p)
 for (i in 1:n){mat[i,]<-as.vector(A[,,i])}
 Xp<-mat%*%(I-(Y1%*%t(Y1)))
 Xp1<-Xp+oo
 array(t(Xp1), dim=c(p, k, n))}
####f4.32
tps2d<-function(M, matr, matt)
{p<-dim(matr)[1]; q<-dim(M)[1]; n1<-p+3
 P<-matrix(NA, p, p)
 for (i in 1:p)
 {for (j in 1:p){
   r2<-sum((matr[i,]-matr[j,])^2)
   P[i,j]<- r2*log(r2)}}
 P[which(is.na(P))]<-0
 Q<-cbind(1, matr)
 L<-rbind(cbind(P,Q), cbind(t(Q),matrix(0,3,3)))
 m2<-rbind(matt, matrix(0, 3, 2))
 coefx<-solve(L)%*%m2[,1]
 coefy<-solve(L)%*%m2[,2]
 fx<-function(matr, M, coef)
 {Xn<-numeric(q)
  for (i in 1:q)
  {Z<-apply((matr-matrix(M[i,],p,2,byrow=T))^2,1,sum)
   Xn[i]<-coef[p+1]+coef[p+2]*M[i,1]+coef[p+3]*M[i,2]+sum(coef[1:p]*(Z*log(Z)))}
  Xn}
 matg<-matrix(NA, q, 2)
 matg[,1]<-fx(matr, M, coefx)
 matg[,2]<-fx(matr, M, coefy)
 matg}
####f4.33
#modified (14.03.14)
tps2<- function(matr, matt, n, master){
  xm<-min(matt[,1])
  ym<-min(matt[,2])
  xM<-max(matt[,1])
  yM<-max(matt[,2])
  rX<-xM-xm; rY<-yM-ym
  a<-seq(xm-1/5*rX, xM+1/5*rX, length=n)
  b<-seq(ym-1/5*rX, yM+1/5*rX,by=(xM-xm)*7/(5*(n-1)))
  m<-round(0.5+(n-1)*(2/5*rX+ yM-ym)/(2/5*rX+ xM-xm))
  M<-as.matrix(expand.grid(a,b))
  ngrid<-tps2d(M,matr,matt)
  r1<-range(master[,1,]); r2<-range(master[,2,])
  r1.1<- abs(r1[2]-r1[1]); r2.1<- abs(r2[2]-r2[1])
  xlim<- c(r1[1]-r1.1/4, r1[2]+r1.1/4)
  ylim<- c(r2[1]-r2.1/4, r2[2]+r2.1/4)
  plot(ngrid, cex=0.2,asp=1,axes=F,xlab="",ylab="", col="gray", xlim=xlim, ylim=ylim)
  for (i in 1:m){lines(ngrid[(1:n)+(i-1)*n,], col="gray")}
  for (i in 1:n){lines(ngrid[(1:m)*n-i+1,], col="gray")}}
####f5.4
fourier1<-function(M, n)
{p<-dim(M)[1]
 an<-numeric(n)
 bn<-numeric(n)
 Z<-complex(real=M[,1],imaginary=M[,2])
 r<-Mod(Z)
 angle<-Arg(Z)
 ao<- 2* sum(r)/p
 for (i in 1:n){
   an[i]<-(2/p)*sum(r * cos(i*angle))
   bn[i]<-(2/p)*sum(r * sin(i*angle))}
 list(ao=ao, an=an, bn=bn )}
####f5.5
ifourier1<-function(ao, an, bn, n, k)
{theta<-seq(0,2*pi, length=n+1)[-(n+1)]
 harm <- matrix (NA, k, n)
 for (i in 1:k)
 {harm[i,]<-an[i]*cos(i*theta)+ bn[i]*sin(i*theta)}
 r<-(ao/2) + apply(harm, 2, sum)
 Z<-complex(modulus=r, argument=theta)
 list(angle = theta, r = r, X = Re(Z), Y = Im(Z))}
####f5.6
fourier2<-function(M,n)
{p<-dim(M)[1]
 an<-numeric(n)
 bn<-numeric(n)
 tangvect<-M-rbind(M[p,],M[-p,])
 perim<-sum(sqrt(apply((tangvect)^2, 1, sum)))
 v0<-(M[1,]-M[p,])
 tet1<-Arg(complex(real=tangvect[,1],imaginary = tangvect [,2]))
 tet0<-tet1[1]
 t1<-(seq(0, 2*pi, length= (p+1)))[1:p]
 phi<-(tet1-tet0-t1)%%(2*pi)
 ao<- 2* sum(phi)/p
 for (i in 1:n){
   an[i]<- (2/p) * sum( phi * cos (i*t1))
   bn[i]<- (2/p) * sum( phi * sin (i*t1))}
 list(ao=ao, an=an, bn=bn, phi=phi, t=t1, perimeter=perim, thetao=tet0)}
####f5.7
ifourier2<-function(ao,an,bn,n,k, thetao=0)
{theta<-seq(0,2*pi, length=n+1)[-(n+1)]
 harm <- matrix (NA, k, n)
 for (i in 1:k)
 {harm[i,]<-an[i]*cos(i*theta)+bn[i]*sin(i*theta)}
 phi<-(ao/2) + apply(harm, 2, sum)
 vect<-matrix(NA,2,n)
 Z<-complex(modulus=(2*pi)/n,argument=phi+theta+thetao)
 Z1<-cumsum(Z)
 list(angle=theta, phi=phi, X=Re(Z1), Y=Im(Z1))}
####f5.8
efourier<-function(M, n=dim(M)[1]/2)
{p<-dim(M)[1]
 Dx<-M[,1]-M[c(p,(1:p-1)),1]
 Dy<-M[,2]-M[c(p,(1:p-1)),2]
 Dt<-sqrt(Dx^2+Dy^2)
 t1<-cumsum(Dt)
 t1m1<-c(0, t1[-p])
 T<-sum(Dt)
 an<-bn<-cn<-dn<-numeric(n)
 for (i in 1:n){
   an[i]<- (T/(2*pi^2*i^2))*sum((Dx/Dt)*(cos(2*i*pi*t1/T)-cos(2*pi*i*t1m1/T)))
   bn[i]<- (T/(2*pi^2*i^2))*sum((Dx/Dt)*(sin(2*i*pi*t1/T)-sin(2*pi*i*t1m1/T)))
   cn[i]<- (T/(2*pi^2*i^2))*sum((Dy/Dt)*(cos(2*i*pi*t1/T)-cos(2*pi*i*t1m1/T)))
   dn[i]<- (T/(2*pi^2*i^2))*sum((Dy/Dt)*(sin(2*i*pi*t1/T)-sin(2*pi*i*t1m1/T)))}
 ao<-2*sum(M[,1]*Dt/T)
 co<-2*sum(M[,2]*Dt/T)
 list(ao=ao,co=co,an=an,bn=bn,cn=cn,dn=dn)}
####f5.9
iefourier<-function(an,bn,cn,dn,k,n,ao=0,co=0)
{theta<-seq(0,2*pi, length=n+1)[-(n+1)]
 harmx <- matrix (NA, k, n)
 harmy <- matrix (NA, k, n)
 for (i in 1:k){
   harmx[i,]<-an[i]*cos(i*theta)+bn[i]*sin(i*theta)
   harmy[i,]<-cn[i]*cos(i*theta)+dn[i]*sin(i*theta)}
 x<-(ao/2) + apply(harmx, 2, sum)
 y<-(co/2) + apply(harmy, 2, sum)
 list(x=x, y=y)}
####f5.10. this function has been sligtly modified (see the errata)
NEF<-function(M, n=dim(M)[1]/2,start=F)
{ef<-efourier(M,n)
 A1<-ef$an[1]; B1<-ef$bn[1]
 C1<-ef$cn[1]; D1<-ef$dn[1]
 theta<-(0.5*atan(2*(A1*B1+C1*D1)/(A1^2+C1^2-B1^2-D1^2)))%%pi
 phaseshift<-matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),2,2)
 M2<-matrix(c(A1,C1,B1,D1),2,2)%*%phaseshift
 v<-apply(M2^2,2, sum)
 if (v[1]<v[2]){theta<-theta+pi/2} 
 theta<-(theta+pi/2)%%pi-pi/2
 Aa<-A1*cos(theta)+B1*sin(theta)
 Cc<-C1*cos(theta)+D1*sin(theta)
 scale<-sqrt(Aa^2+Cc^2)
 psi<-atan(Cc/Aa)
 if (Aa<0){psi<-psi+pi}
 size<-(1/scale)
 rotation<-matrix(c(cos(psi),-sin(psi),sin(psi),cos(psi)),2,2)
 A<-B<-C<-D<-numeric(n)
 if (start){theta<-0}
 for (i in 1:n){
   mat<-size*rotation%*%matrix(c(ef$an[i],ef$cn[i],ef$bn[i],ef$dn[i]),2,2)%*%matrix(c(cos(i*theta),sin(i*theta),-sin(i*theta),cos(i*theta)),2,2)
   A[i]<-mat[1,1]
   B[i]<-mat[1,2]
   C[i]<-mat[2,1]
   D[i]<-mat[2,2]}
 list(A=A,B=B,C=C,D=D,size=scale,theta=theta,psi=psi,ao=ef$ao,co=ef$co, Aa=Aa, Cc=Cc)}

#-------------------------------------------------------------
# Generate folds
# @description Give the folds of k-fold cross-validation
# @author Matthias C. M. Troffaes
# @references https://gist.github.com/mcmtroffaes/709908
#-------------------------------------------------------------

kfcv.sizes = function(n, k=10) {
  sizes = c()
  for (i in 1:k) {
    first = 1 + (((i - 1) * n) %/% k)
    last = ((i * n) %/% k)
    sizes = append(sizes, last - first + 1)
  }
  sizes
}

kfcv.testing = function(n, k=10) {
  indices = list()
  sizes = kfcv.sizes(n, k=k)
  values = 1:n
  for (i in 1:k) {
    s = sample(values, sizes[i])
    indices[[i]] = s
    values = setdiff(values, s)
  }
  indices
}

#-------------------------------------------------------------
# aligne configuration
# @description aligne the landmark configuration to the PC axis
# @details modified from J Claude's Book Function, which in the book, 
#  the \code{aligne} function takes array, here, the \code{aligne2} function
# takes matrix
# @references Claude J. (2008). Morphometrics with R. Springer
#-------------------------------------------------------------

aligne2<-function(A){
  B<-A
  Ms<-scale(A, scale=F)
  sv<-eigen(var(Ms))
  M<-Ms%*%sv$vectors
  B<-M
  B
} 