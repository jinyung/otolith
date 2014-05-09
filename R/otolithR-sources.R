# ====== FUNCTIONS SOURCED FROM OTHER PLACES ====== #
# 1. JClaude Book functions                   (line 12)
# 2. kfcv.sizes() and kfcv.testing()          (line 903) 
# 3. aligne2()                                (line 932)
#
#
#
#
#

#-------------------------------------------------------------
# JClaude Book functions
# @author Julien Claude
# @references Claude J. (2008). Morphometrics with R. Springer
# @description functions developped in the book.
# @note \code{tps} function is modified
#-------------------------------------------------------------

######functions developped in Claude J. (2008). Morphometrics with R. Springer. 
######Refer to the book if you use them
####f1.1
myfun<-function(vec)
{ sum(vec*vec) }
####f2.1
Conte<-function(x, imagematrix)
{I<-imagematrix
 x<-rev(x)
 x[1]<-dim(I)[1]-x[1]
 while (abs(I[x[1],x[2]]-I[x[1],(x[2]-1)])<0.1){x[2]<-x[2]-1}
 a<-1
 M<-matrix(c(0,-1,-1,-1,0,1,1,1,1,1,0,-1,-1,-1,0,1),2,8,byrow=T)
 M<-cbind(M[,8],M,M[,1])
 X<-0; Y<-0;
 x1<-x[1]; x2<-x[2]
 SS<-NA; S<-6
 while ((any(c(X[a],Y[a])!=c(x1,x2) ) | length(X)<3))
 {if (abs(I[x[1]+M[1,S+1],x[2]+M[2,S+1]]-I[x[1],x[2]])<0.1)
 {a<-a+1;X[a]<-x[1];Y[a]<-x[2];x<-x+M[,S+1]
  SS[a]<-S+1; S<-(S+7)%%8}
 else if (abs(I[x[1]+M[1,S+2],x[2]+M[2,S+2]]-I[x[1],x[2]])<0.1)
 {a<-a+1;X[a]<-x[1];Y[a]<-x[2];x<-x+M[,S+2]
  SS[a]<-S+2; S<-(S+7)%%8}
 else if (abs(I[x[1]+M[1,(S+3)],x[2]+M[2,(S+3)]]-I[x[1],x[2]])<0.1)
 {a<-a+1;X[a]<-x[1];Y[a]<-x[2];x<-x+M[,(S+3)]
  SS[a]<-S+3; S<-(S+7)%%8}
 else S<-(S+1)%%8}
 list(X=(Y[-1]), Y=((dim(I)[1]-X))[-1])}
####f2.2
ild<-function(E, F){ sqrt(sum((E-F)^2))}
####f2.3
angle<-function(v1, v2)
{temp <- sum(v1*v2)/( sqrt(sum(v1^2))*sqrt(sum(v2^2)) )
 acos (temp)}
####f2.4
angle3<-function(v1, v2)
{a<-angle(v1, v2)
 b<-sign( det(rbind(1, v1, v2)) )
 if (a == 0 & b == 1){jo<-pi/2}
 else if (a == 0 & b == -1){jo<- - pi/2}
 else {jo<- a * b}
 (pi+jo)%%(2*pi)-pi}
####f2.5
regularradius<-function(Rx, Ry, n)
{le<-length(Rx)
 M<-matrix(c(Rx, Ry), le, 2)
 M1<-matrix(c(Rx-mean(Rx), Ry-mean(Ry)), le, 2)
 V1<-complex(real=M1[,1], imaginary=M1[,2])
 M2<-matrix(c(Arg(V1), Mod(V1)), le, 2)
 V2<-NA
 for (i in 0:(n-1))
 {V2[i+1]<-which.max((cos(M2[,1]-2*i*pi/n)))}
 V2<-sort(V2)
 list("pixindices"=V2,"radii"=M2[V2,2],"coord"=M1[V2,])}
####f2.6
smoothout<-function(M, n)
{p<-dim(M)[1]
 a<-0
 while (a<=n)
 {a<-a+1
  Ms<-rbind(M[p,],M[-p,])
  Mi<-rbind(M[-1,],M[1,])
  M<-M/2+Ms/4+Mi/4}
 M}
####f2.7
landmark.addition<-function(M, n)
{a<-0
 while(a<=n)
 {p<-dim(M)[1]
  k<-dim(M)[2]
  N<-matrix (NA,2*p,k)
  N[((1:p)*2)-1,]<-M
  N[(1:p)*2,]<-(M+(rbind(M[-1,],M[1,])))/2
  M<-N}
 M}
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
####f3.1
ELLI<-function(x,y,conf=0.95,np=50)
{centroid<-apply(cbind(x,y),2,mean)
 ang <- seq(0,2*pi,length=np)
 z<-cbind(cos(ang),sin(ang))
 radiuscoef<-qnorm((1-conf)/2, lower.tail=F)
 vcvxy<-var(cbind(x,y))
 r<-cor(x,y)
 M1<-matrix(c(1,1,-1,1),2,2)
 M2<-matrix(c(var(x), var(y)),2,2)
 M3<-matrix(c(1+r, 1-r),2,2, byrow=T)
 ellpar<-M1*sqrt(M2*M3/2)
 t(centroid + radiuscoef * ellpar %*% t(z))}
####f3.2#modif line 106, 109, 114, 116 #angle varie from 0 to 90°
coli<-function(ev1, ev2, nperm=1000, graph=T)
{dist<-numeric(nperm)
 n<-length(ev1)
 Angle<-function(v1, v2)
 {90-abs((180*(acos(sum(v1*v2)/(sqrt(sum(v1^2))*sqrt(sum(v2^2)))))/pi)-90)}
 for (i in 1:nperm)
 {X1<-runif(n, -1, 1); X2<-runif(n, -1, 1)
  dist[i]<-Angle(X1, X2)}
 zobs<-Angle(ev1, ev2)
 pv<-length(dist[dist<zobs])/nperm
 if (graph)
 {hist(dist,breaks=50,
       main="Distribution of angles between 2 random vectors", xlab="Z statistic",ylab="# of random vect.", sub=paste("Actual z-obs =",round(zobs,5),": p<",round(pv),5))
  abline(v=zobs)}
 list(angle=zobs,p=pv)}
####f3.3
isojoli<-function(mat)
{n<-dim(mat)[1]; p<-dim(mat)[2]
 S<-var(log(mat))
 V1<-rep(sqrt(1/p),p)
 L1<-svd(S)$d[1]
 chiobs<-(n-1)*(L1*t(V1)%*%ginv(S)%*%V1+(1/L1)*t(V1)%*%S%*%V1-2)
 unlist(list(Chisq=chiobs, p=pchisq(chiobs,p-1,lower.tail=F)))}
####f3.4
pls<-function(M1, M2)
{p1<-dim(M1)[2]; p2<-dim(M2)[2]; n<-dim(M1)[1]
 sM12<-svd(var(cbind(M1,M2))[1:p1, (p1+1):(p1+p2)])
 vM12<-var(cbind(M1,M2))[1:p1, (p1+1):(p1+p2)]
 vM21<-var(cbind(M1,M2))[(p1+1):(p1+p2), 1:p1]
 v11<-var(M1)
 v22<-var(M2)
 D<-sM12$d; F1<-sM12$u; F2<-sM12$v
 Rv<-sum(diag(vM12%*%vM21))/sqrt(sum(diag(v11%*%v11))*sum(diag(v22%*%v22)))
 list(Rv=Rv, F1=F1, F2=F2, D=D)}
####f4.1
centcoord<-function(M){apply(M,2,mean)}
####f4.2
centsiz<-function(M)
{p<-dim(M)[1]
 size<-sqrt(sum(apply(M, 2,var))*(p-1))
 list("centroid_size" = size,"scaled" = M/size)}
####f4.3
basesiz<-function(M, p1, p2)
{sqrt(sum((M[p1,]-M[p2,])^2))}
####f4.4
booksteinM<-function(M, p1, p2)
{D<-basesiz(M, p1, p2)
 m<-matrix(NA, nrow(M), ncol(M))
 p1<-M[p1,]
 p2<-M[p2,]
 m[,1]<-(((p2[1]-p1[1])*(M[,1]-p1[1])+(p2[2]-p1[2])*(M[,2]-p1[2]))/(D^2))-0.5
 m[,2]<-((p2[1]-p1[1])*(M[,2]-p1[2])-(p2[2]-p1[2])*(M[,1]-p1[1]))/(D^2)
 m}
####f4.5
booksteinA<-function(A, p1, p2)
{B<-array(NA, dim=c(dim(A)[1],dim(A)[2],dim(A)[3]))
 for (i in 1: dim(A)[3])
 {B[,,i]<-booksteinM(A[,,i], p1, p2)}
 B}
####f4.6
mbshape<-function(A,p1, p2)
{B<-booksteinA(A, p1, p2)
 k<-dim(A)[2]
 mbshape<-matrix(NA, dim(A)[1], dim(A)[2])
 for (i in 1:k)
 {mbshape[,i]<-apply(B[,i,], 1, mean)}
 mbshape}
####f4.7
tranb<-function(M, p1, p2)
{M-matrix((M[p1,]+M[p2,])/2, nrow(M), ncol(M), byrow=T)}
####f4.8
transb<-function(M, p1, p2)
{tranb(M, p1, p2)/basesiz(M, p1, p2)}
####f4.9
bookstein3d<-function(M, p1, p2, p3)
{m<-transb(M, p1, p2)
 te1<-angle3(c(1,0,0),c(m[p2,1],m[p2,2],0))
 Rz<-matrix(c(cos(te1),sin(te1),0,-sin(te1),cos(te1),0,0,0,1),3,3,byrow=T)
 m1<-t(Rz%*%t(m))
 om1<-angle3(c(m1[p2,1],0,0),c(m1[p2,1],0,m1[p2,3]))
 Ry<-matrix(c(cos(om1),0,-sin(om1),0,1,0,sin(om1),0,cos(om1)),3,3,byrow=T)
 m2<-t((Ry)%*%(Rz)%*%t(m))/basesiz(M,p1,p2)
 ph1<-angle3(c(0,m2[p3,2],0),c(0,m2[p3,2:3]))
 Rx<-matrix(c(1,0,0,0,cos(ph1),sin(ph1),0,-sin(ph1),cos(ph1)),3,3,byrow=T)
 t(Rx%*%Ry%*%Rz%*%t(m))}
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
####f4.23
stp<-function(A)
{p<-dim(A)[1];k<-dim(A)[2];n<-dim(A)[3]
 Yn<-mshape(A)
 B<-array(NA, dim=c(p,k,n))
 for (i in 1:n)
 {rho<-2*asin((sqrt(sum((A[,,i]-Yn)^2)))/2)
  B[,,i]<-A[,,i]/(cos(rho))}
 return(B)}
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
####f4.25
medsize<-function(M)
{mat<-as.matrix(dist(M))
 median(apply(mat, 2, median, na.rm=TRUE))}
####f4.26
argallvec<-function(X1, X2)
{p<-dim(X1)[1]
 m<-m<-matrix(NA, p, p)
 for (i in 1:p){ for (j in 1:p){
   m[i,j]<-Arg(complex(1,X2[i,1],X2[i,2])-complex(1,X2[j,1],X2[j,2]))-Arg(complex(1,X1[i,1],X1[i,2])-complex(1,X1[j,1],X1[j,2]))  }}
 ((m+pi)%%(2*pi))-pi}
####f4.27
oPsup<-function(M1, M2)
{p<-dim(M1)[1]; k<-dim(M1)[2]
 ols<-fPsup(M1, M2)
 X1<-ols[[1]]; X2<-ols[[2]]
 M<-as.matrix(dist(X2))/as.matrix(dist(X1))
 beta<-median(apply(M, 2, median, na.rm=TRUE))
 ARG<-argallvec(X1,X2)
 phi<-median(apply(ARG, 2,median, na.rm=TRUE))
 Gam<-matrix(c(cos(phi),-sin(phi),sin(phi),cos(phi)),2,2)
 alpha<-X2-beta*X1%*%Gam
 alpha<-matrix(apply(alpha,2,median), p, k, byrow=T)
 list("Mo1"=beta*X1%*%Gam + alpha,"Mo2"=X2)}
####f4.28
rrot<-function(M1, M2)
{p<-dim(M1)[1]; k<-dim(M1)[2]
 P<-1:p
 I<-diag(1,3)
 A<-array(NA, dim=c(p,p,p,3,3))
 for (i in P){
   for (j in P[-i]){
     for (k in P[-c(i,j)]){
       M1i<-matrix(M1[i,])
       M1j<-matrix(M1[j,])
       M1k<-matrix(M1[k,])
       M2i<-matrix(M2[i,])
       M2j<-matrix(M2[j,])
       M2k<-matrix(M2[k,])
       M1ij<-(M1j-M1i)/sqrt(sum((M1j-M1i)^2))
       M2ij<-(M2j-M2i)/sqrt(sum((M2j-M2i)^2))
       M1ijk<-((M1k-M1i)-as.numeric(t(M1k-M1i)%*%M1ij)*M1ij)/ sqrt(sum(((M1k-M1i)-as.numeric(t(M1k-M1i)%*%M1ij)*M1ij)^2))
       M2ijk<-((M2k-M2i)-as.numeric(t(M2k-M2i)%*%M2ij)*M2ij)/ sqrt(sum(((M2k-M2i)-as.numeric(t(M2k-M2i)%*%M2ij)*M2ij)^2))
       M1s<-matrix(c(M1ij[2]*M1ijk[3]-M1ij[3]*M1ijk[2],M1ij[3]*M1ijk[1]-M1ij[1]*M1ijk[3],M1ij[1]*M1ijk[2]-M1ij[2]*M1ijk[1]))
       M2s<-matrix(c(M2ij[2]*M2ijk[3]-M2ij[3]*M2ijk[2],M2ij[3]*M2ijk[1]-M2ij[1]*M2ijk[3],M2ij[1]*M2ijk[2]-M2ij[2]*M2ijk[1]))
       X<-cbind(M1ij, M1ijk, M1s)
       Y<-cbind(M2ij, M2ijk, M2s)
       H<-t(X)%*%Y
       A[i,j,k,,]<-S<-solve(H+I)%*%(H-I)}}}
 S1<-apply(apply(apply(A,c(1,2,4,5),median,na.rm=T),c(1,3,4),median,na.rm=T),2:3,median,na.rm=T)
 H<-(I+S1)%*%solve(I-S1)
 list("ref"=M2, "targ"=M1%*%t(H), "Gamma"=H)}
####f4.29
r3sup<-function(M1, M2)
{p<-dim(M1)[1]; k<-dim(M1)[2]
 ols<-fPsup(M1, M2)
 X1<-ols[[1]]; X2<-ols[[2]]
 B<-as.matrix(dist(X2))/as.matrix(dist(X1))
 M<-B; M<-t(M)
 M[row(B)<col(B)]<-B[row(B)<col(B)]
 M<-as.matrix(dist(X2))/as.matrix(dist(X1))
 beta<-median(apply(M, 2, median, na.rm=TRUE))
 Gam<-rrot(X1,X2)[[3]]
 alpha<-X2-beta*X1%*%Gam
 alpha<-matrix(apply(alpha,2,median),p,k, byrow=T)
 list("Mo1"=beta*X1%*%Gam + alpha,"Mo2"=X2)}
####f4.30
grf2<-function(A)
{p<-dim(A)[1]; k<-dim(A)[2]; n<-dim(A)[3]
 A<-pgpa(A)$rotated
 D<-B<-array(NA, dim=c(p,k,n))
 for (i in 1:n){
   B[,,i]<-A[,,i]/medsize(A[,,i])}
 Y<-apply(B, 1:2, median)
 Y<-Y/medsize(Y)
 A0<-10
 iter<-1
 while(A0>0.0005){
   for (i in 1:n){
     M<-as.matrix(dist(Y))/as.matrix(dist(B[,,i]))
     beta<-median(apply(M, 2, median, na.rm=TRUE))
     ARG<-argallvec(B[,,i],Y)
     phi<-median(apply(ARG, 2,median, na.rm=TRUE))
     Gam<-matrix(c(cos(phi),-sin(phi),sin(phi),cos(phi)),2,2)
     alpha<-Y-beta*B[,,i]%*%Gam
     D[,,i]<-beta*B[,,i]%*%Gam+matrix(apply(alpha,2,median),p,k,byrow=T)}
   Yb<-apply(D, 1:2, median)
   A0<-median(sqrt(apply((Yb-Y)^2,1,sum)))
   Y<-Yb<-Yb/medsize(Yb)
   B<-D
   iter<-iter+1}
 list("rotated"=D,"limit"=A0,"iteration"=iter,"medshape"=Yb)}
####f4.31
grf3<-function(A)
{p<-dim(A)[1]; k<-dim(A)[2]; n<-dim(A)[3]
 A<-pgpa(A)$rotated
 D<-B<-array(NA, dim=c(p,k,n))
 for (i in 1:n){
   B[,,i]<-A[,,i]/medsize(A[,,i])}
 Y<-apply(B, 1:2, median)
 Y<-Y/medsize(Y)
 A0<-10
 iter<-1
 while(A0>0.005){
   for (i in 1:n){
     M<-as.matrix(dist(Y))/as.matrix(dist(B[,,i]))
     beta<-median(apply(M, 2, median, na.rm=TRUE))
     Gam<-rrot(B[,,i],Y)[[3]]
     alpha<-Y-beta*B[,,i]%*%Gam
     D[,,i]<-beta*B[,,i]%*%Gam + matrix(apply(alpha,2,median),p,k,byrow=T)}
   Yb<-apply(D, 1:2, median)
   A0<-median(sqrt(apply((Yb-Y)^2,1,sum)))
   Y<-Yb<-Yb/medsize(Yb)
   B<-D
   iter<-iter+1}
 list("rotated"=D,"limit"=A0,"iteration"=iter,"medshape"=Yb)}
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


####f4.34
FM<-function(M){as.matrix(dist(M))}
####f4.34'
fm<-function(M){mat<-FM(M); mat[col(mat)<row(mat)]}
####f4.35
mEDMA<-function(A)
{n<-dim(A)[3];p<-dim(A)[1]; k<-dim(A)[2]
 E<-matrix(NA,n,p*(p-1)/2)
 for (i in 1:n){E[i,]<-(fm(A[,,i]))^2}
 Em<-apply(E,2,mean)
 S<-(apply(t((t(E)-Em)^2),2,sum))/n
 if (k==2){omega<-(Em^2-S)^0.25}
 if (k==3){omega<-(Em^2-1.5*S)^0.25}
 Om<-diag(0,p)
 Om[row(Om)>col(Om)]<-omega; Om<-t(Om)
 Om[row(Om)>col(Om)]<-omega
 Om}
####f4.36
MDS<-function(mat, k)
{p<-dim(mat)[1]
 C1<-diag(p)-1/p*matrix(1,p,p)
 B<- -0.5*C1%*%mat^2%*%C1
 eC<-eigen(B)
 eve<-eC$vectors
 eva<-eC$values
 MD<-matrix(NA, p, k)
 for (i in 1:k)
 {MD[,i]<-sqrt(eva[i])*eve[,i]}
 MD}
####f4.37
mEDMA2<-function(A)
{k<-dim(A)[2]
 Eu<-mEDMA(A)
 M<-MDS(Eu,k)
 list("M"=M, "FM"=FM(M))}
####f4.38
vEDMA<-function(A)
{p<-dim(A)[1]
 k<-dim(A)[2]; n<-dim(A)[3]
 Bs<-array(NA, dim=c(p,p,n))
 for (i in 1:n){
   Cc<-apply(A[,,i],2,mean)
   Ac<-t(t(A[,,i])-Cc)
   Bs[,,i]<-Ac%*%t(Ac)}
 B<-apply(Bs, 1:2,mean)
 M<-mEDMA2(A)$M
 Ek<-(B-M%*%t(M))/k
 Ek}
####f4.39
alltri<-function(p)
{t1<-expand.grid(1:p,1:p,1:p)
 t2<-t1[t1[,1]<t1[,2],]
 t3<-t2[t2[,2]<t2[,3],]
 t3}
####f4.40
anglerao<-function(M, triang)
{n<-dim(triang)[1]
 triangl<-triang
 for (i in 1:n){
   t1<-triang[i,]
   AB<-complex(real=M[t1[2],1]-M[t1[1],1],imag=M[t1[2],2]-M[t1[1],2])
   BC<-complex(real=M[t1[3],1]-M[t1[2],1],imag=M[t1[3],2]-M[t1[2],2])
   AC<-complex(real=M[t1[3],1]-M[t1[1],1],imag=M[t1[3],2]-M[t1[1],2])
   BAC<-Arg(complex(argument=-Arg(AB)+Arg(AC)))
   CBA<-Arg(complex(argument=-Arg(BC)+Arg(-AB)))
   ACB<-Arg(complex(argument=-Arg(-AC)+Arg(-BC)))
   triangl[i,]<-c(BAC,CBA,ACB)}
 triangl}
####f4.41
meanrao<-function(A,triang)
{n<-dim(A)[3];p<-dim(A)[1];nt<-dim(triang)[1]
 A1<-array(NA, dim=c(nt, 3, n))
 A2<-matrix(NA, nt, 3)
 A3<-matrix(0, nt, 3)
 for(i in 1:n)
 {A1[,,i]<-anglerao(A[,,i], triang)}
 for(j in 1:nt)
 {Aj<-A1[j,,]
  mAj<-apply(Aj,1,mean)
  MAj<-c(max(abs(Aj[1,])),max(abs(Aj[2,])),max(abs(Aj[3,])))
  MAji<-which(MAj==max(MAj))
  A3[j,MAji]<-1
  mAj[MAji]<-(pi-sum(abs(mAj[-MAji])))*sign(sum(mAj[-MAji]))
  A2[j,]<-mAj}
 list(mean=A2, maxangle=A3)}
####f4.42
vcvrao<-function(A,triang)
{n<-dim(A)[3];nt<-dim(triang)[1]
 Aa<-t(matrix(A, nt*3,n))
 VCV<-var(Aa)
 maxangle<-meanrao(A,triang)$maxangle
 strip<-which(maxangle==1)
 VCV<-VCV[-strip, -strip]
 list(VCV=VCV,maxangle=maxangle, strip=strip)}
####f4.43
raoinv<-function(T, triang)
{nt<-dim(triang)[1]
 M<-NA
 triang1<-triang
 M1<-M[triang[1,1]]<-0+0i
 M2<-M[triang[1,2]]<-1+0i
 M[triang[1,3]]<-complex(modulus=sin(T[1,2])*Mod(M2-M1)/sin(T[1,3]),argument=T[1,1]+Arg(M2-M1))+M[triang[1,1]]
 itriang<-triang1[1,]
 triang1<-triang1[-1,]
 T1<-T[-1,]
 i<-1
 while (length(itriang)<length(unique(as.vector(triang))))
 {while (length(which(itriang%in%triang1[i,]==TRUE))<2)
 {i<-i+1}
 iP<-which(itriang%in%triang1[i,]==TRUE)
 P<-itriang[iP]
 iN<-which(triang1[i,]%in%itriang==FALSE)
 N<-triang1[i,iN]
 P1<-P[1]; iP1<-which(triang1[i,]==P1)
 P2<-P[2]; iP2<-which(triang1[i,]==P2)
 if (iP1==1 & iP2==3){P3<-P1; P1<-P2; P2<-P3; iP1<-3; iP2<-1}
 if (iP1==2 & iP2==1){P3<-P1; P1<-P2; P2<-P3; iP1<-1; iP2<-2}
 if (iP1==3 & iP2==2){P3<-P1; P1<-P2; P2<-P3; iP1<-2; iP2<-3}
 M1<-M[P1]
 M2<-M[P2]
 M[N]<-complex(modulus=sin(T1[i,iP2])*Mod(M2-M1)/sin(T1[i,iN]),argument=T1[i,iP1]+Arg(M2-M1))+M[P1]
 T1<-T1[-i,]
 itriang<-unique(c(itriang, triang1[i,]))
 triang1<-triang1[-i,]
 i<-1}
 cbind(Re(M), Im(M))}
####f5.1
cumchord<-function(M)
{cumsum(sqrt(apply((M-rbind(M[1,],M[-(dim(M)[1]),]))^2,1,sum)))}
####f5.2
bezier<-function(M,n = dim(M)[1])
{p<-dim(M)[1]
 if (n != p) {n<-n+1}
 M1<-M/cumchord(M)[p]
 t1<-1-cumchord(M1)
 J<-matrix(NA, p, p)
 for (i in 1:p){
   for (j in 1:p){J[i, j]<-(factorial(p-1)/(factorial(j-1)*factorial(p-j)))*(((1-t1[i])^(j-1))*t1[i]^(p-j))}}
 B<-ginv(t(J[,1:n])%*%J[,1:n])%*%(t(J[,1:n]))%*%M
 M<-J[,1:n]%*%B
 list(J=J, B=ginv(t(J[,1:n])%*%J[,1:n])%*%(t(J[,1:n]))%*%M)}
####f5.3
beziercurve<-function(B,p)
{X<-Y<-numeric(p)
 n<-dim(B)[1]-1
 t1<-seq(0,1,length=p)
 coef<-choose(n, k=0:n)
 b1<-0:n; b2<-n:0
 for (j in 1:p)
 {vectx<-vecty<-NA
  for (i in 1:(n+1))
  {vectx[i]<-B[i,1]*coef[i]*t1[j]^b1[i]*(1-t1[j])^b2[i]
   vecty[i]<-B[i,2]*coef[i]*t1[j]^b1[i]*(1-t1[j])^b2[i]}
  X[j]<-sum(vectx)
  Y[j]<-sum(vecty)}
 cbind(X, Y)}
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
####f5.11
eigenshape<-function(Z)
{n<-dim(Z)[1]
 Zz<-scale(Z, scale=F)
 R<-t(Zz)%*%Zz
 V<-svd(R)$u
 S<-svd(R)$d
 U<-Zz%*%V%*%(diag(1/sqrt(S)))
 weight<-diag(sqrt(S))%*%t(V)
 list(eigenshape=U, weights=weight[1:n,], eigenvalues=S)}
####f6.1
procalign<-function(A)
{pA<-pgpa(A); n<-dim(A)[3]; k<-dim(A)[2]
 A<-pA$rotated; msh<-pA$mshape
 A1<-A
 sv<-eigen(var(msh))
 V<-sv$vectors
 rotmsh<-msh%*%V
 for (i in 1:n){A1[,,i]<-pPsup(A[,,i],rotmsh)$Mp1}
 list("rotated"=orp(A1), "meansh"=rotmsh)}
####f6.2
uniform2D<-function(A)
{n<-dim(A)[3]
 kp<-dim(A)[1]*dim(A)[2]
 temp<-procalign(A)
 msh<-temp$meansh
 proc<-temp$rotated
 X<-t(matrix(proc,kp,n))
 V<-X-rep(1, n)%*%t(as.vector(msh))
 alph<-sum(msh[,1]^2)
 gam<-sum(msh[,2]^2)
 U1<-c(sqrt(alph/gam)*msh[,2],sqrt(gam/alph)*msh[,1])
 U2<-c(-sqrt(gam/alph)*msh[,1],sqrt(alph/gam)*msh[,2])
 score<-V%*%cbind(U1,U2)
 list("scores"=score,"uniform"=cbind(U1,U2),"meanshape"=msh,"rotated"=proc)}
####f6.3
uniformG<-function(A)
{n<-dim(A)[3]; k<-dim(A)[2]; p<-dim(A)[1]
 temp<-procalign(A)
 msh<-temp$meansh
 proc<-temp$rotated
 Bn<-Xn<-list()
 for (i in 1:k)
 {Bn[[i]]<-solve(t(msh)%*%msh)%*%t(msh)%*%(proc[,i,])}
 if (k==2)
 {LSR<-svd(t(rbind(Bn[[1]],Bn[[2]]))%*%(diag(1,2)%x%t(msh)))
  score<-LSR$u%*%diag(LSR$d)[,2:3]
  Un<-LSR$v[,2:3]}
 if (k==3)
 {LSR<-svd(t(rbind(Bn[[1]],Bn[[2]],Bn[[3]]))%*%(diag(1,3)%x%t(msh)))
  score<-LSR$u%*%diag(LSR$d)[,2:6]
  Un<-LSR$v[,2:6]}
 list("score"=score, "uniform"=Un)}
####f6.4
Hotellingsp<-function(SSef, SSer, dfef, dfer, exact=F)
{library(MASS)
 p <- qr(SSef+SSer)$rank
 k<-dfef; w<-dfer
 s<-min(k,p)
 m<-(w-p-1)/2
 t1<-(abs(p-k)-1)/2
 Ht<-sum(diag(SSef%*%ginv(SSer)))
 Fapprox<-Ht*(2 * (s*m+1))/(s^2*(2*t1+s+1))
 ddfnum<-s*(2*t1+s+1)
 ddfden<-2*(s*m+1)
 pval= 1-pf(Fapprox, ddfnum, ddfden)
 if (exact)
 {b<-(p+2*m)*(k+2*m)/((2*m+1)*(2*m-2))
  c1<-(2+(p*k+2)/(b-1))/(2*m)
  Fapprox<-((4+(p*k+2)/(b-1))/(p*k))*(Ht/c1)
  ddfnum<-p*k
  ddfden<-4+(p*k+2)/(b-1)}
 unlist(list("dfeffect"=dfef,"dferror"=dfer,"T2"=Ht,"Approx_F"=Fapprox,"df1"=ddfnum,"df2"=ddfden,"p"=pval))}
####f6.5
mantrstat<-function (m1, m2)
{ cor(lower.triang(m1),lower.triang(m2))}
####f6.6
permrowscols<-function (m1,n,coord)
{s <- sample(1:(n/coord))
 m1[rep(s,coord), rep(s,coord)]}
####f6.7
mantel.t2<-function(m1,m2,coord=1,nperm=1000,graph=FALSE,...)
{n<-nrow(m1)
 realz<-mantrstat(m1, m2)
 nullstats<-replicate(nperm,mantrstat(m1,permrowscols(m2,n,coord)))
 pval <- sum(nullstats > realz)/nperm
 if (graph) {
   plot(density(nullstats), type = "l", ...)
   abline(v = realz)   }
 list(r.stat = realz, p = pval)}
####f6.8
Rv<-function(VCV, index1)
{VCV1<-VCV[index1,index1]
 VCV2<-VCV[-index1,-index1]
 VCV12<-VCV[index1,-index1]
 VCV21<-VCV[-index1,index1]
 sum(diag(VCV12%*%VCV21))/sqrt(sum(diag(VCV1%*%VCV1))*sum(diag(VCV2%*%VCV2)))}
####f7.1'
lowertriang<-function (m)
{d<- dim(m)
 m[col(m) <= row(m)]}
####f7.1''
mstat<-function (m1, m2)
{cor(lower.triang(m1), lower.triang(m2))}
####f7.1
mantel.t3<-function(m1,m2,coord=1,nperm=5000,graph=FALSE)
{n <- nrow(m1)
 realz <- mstat(m1, m2)
 nullstats<-replicate(nperm,mstat(m1,permrowscols(m2, n, coord)))
 list(r.stat = realz, null=nullstats)}

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