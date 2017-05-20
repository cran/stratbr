optimal_allocation<-function(Nh_Vh,n,LL,minimumsample,X)
{
  Nh=Nh_Vh[1:LL]
  Vh=Nh_Vh[((LL+1):(2*LL))]
  funcao_objetivo=as.vector(outer(1/c(1:n),(Vh*Nh^2)))
  inf=1
  sup=n
  for(i in 1:LL) {funcao_objetivo[inf:sup]=funcao_objetivo[inf:sup]-(Vh[i]*Nh[i]);inf=sup+1;sup=sup+n}
  funcao_objetivo[which(funcao_objetivo<0)]=max(Vh)*max(Nh^2)*10
  restricoes=c(rep("==",LL+1),rep("<=",LL),rep(">=",LL))

  b=c(rep(1,LL),n,Nh,rep(minimumsample,LL))
  A=matrix(0,nrow=(3*LL)+1,ncol=n*LL)
  for(i in 1:LL) {A[i,(1+(n*(i-1))):(i*n)]=rep(1,n)}
  A[LL+1,]=rep(c(1:n),LL)
  inf=1
  sup=n
  for(i in (LL+2):((2*LL)+1)){A[i,inf:sup]=c(1:n);inf=sup+1;sup=sup+n}
  for(i in ((2*LL)+2):((3*LL)+1)) {A[i,]=A[i-LL,]}
  solucao=Rglpk_solve_LP(funcao_objetivo,A,restricoes,b,types=rep("B",length(funcao_objetivo)),max=FALSE)
  status=solucao$status
  var_tot=solucao$optimum
  solucao=solucao$solution
  x=matrix(rep(0,LL),nrow=1,byrow=TRUE)
  x[1,1]=sum(solucao[1:n]*c(1:n))
  if (LL>1)
  {inf=n+1
  sup=2*n
  for(i in 2:LL)
  {x[1,i]=sum(solucao[inf:sup]*c(1:n))
  inf=sup+1
  sup=sup+n
  }
  }
  if (status>0) {var_tot=Inf}
  solution=matrix(c(x,sqrt(var_tot)/X),nrow=1)
  return(solution)
}