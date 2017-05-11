#'Optimization Algorithm to solve stratification problem
#'
#'@param X  Stratification variable.
#'@param stratas  Number of strata.
#'@param nsample  Sample size.
#'@param tampop   Number of chromosomes BRKGA.
#'@param totgen   Number of generations BRKGA.
#'@param pelite   Percentage elite solutions BRKGA.
#'@param pmutant  Percentage mutant solutions BRKGA.
#'@param rc       Crossover probability  BRKGA.
#'@param minimumsample Minimum sample size (smallest possible sample size in any stratum).
#'@param cores Numerical amount of CPUs requested for the cluster.
#'@return \item{cvtot}{Coefficient of variation for the estimator of total of the stratification variable considered.}
#'@return \item{nh}{Number of sample elements, or sample size, in stratum h.}
#'@return \item{Nh}{Number of population elements, or population size, in stratum h.}
#'@return \item{Sh2}{Population variance of the stratification variable x in stratum h.}
#'@return \item{bk}{Strata boundaries}
#'@return \item{cputime}{Time consumed by the algorithm in seconds.}
#'
#'@references Brito, J.A.M, Silva, P.L.N.,Semaan, G.S. and Maculan, N. (2015).
#'            Integer Programming Formulations Applied to Optimal Allocation in Stratified Sampling.
#'            Survey Methodology, 41: 427-442.
#'
#'            Brito, J.A.M, Semaan, G.S., Fadel, A.C. and Brito, L.R.(2017).
#'            An optimization approach applied to the optimal stratification problem,
#'            Communications in Statistics - Simulation and Computation.
#'
#'            Gonçalves, J.R. and Resende, M.G.C. (2011).
#'            Biased random-key genetic algorithms for combinatorial optimization,
#'            Journal of Heuristics, 17: 487-525.
#'
#'@author Jose Brito (jambrito@gmail.com), Pedro Luis and Tomas Veiga.
#'
#'@description This function aims at constructing optimal strata with an optimization algorithm
#'based on a global optimisation technique called Biased Random Key Genetic Algorithms(BRKGA).
#'The optimization algorithm is applied to solve the onedimensional case,
#'which reduces the stratification problem to just determining strata boundaries.
#'Assuming that the number H of strata and the total sample size n are fixed,
#'it is possible to produce the strata boundaries by taking into consideration
#'an objective function associated with the variance.
#'This function determines strata boundaries so that the elements in each
#'stratum are more homogeneous among themselves.
#'@export
#'@import Rglpk
#'@import snowfall
#'@import stratification
#'@importFrom stats runif
#'@importFrom stats var
#'@importFrom stats quantile
#'@importFrom utils flush.console
#'@examples
#'data(Sweden)
#'REV84<-Sweden[,9]
#'solution1<-stratbr(REV84,stratas=3,nsample=100,cores=4,totgen=20)
#'#'solution2<-stratbr(REV84,stratas=3,nsample=100,cores=4,totgen=20,minimumsample=30)
#'ME84<-Sweden[,8]
#'solution3<-stratbr(ME84,stratas=4,nsample=75,cores=2,totgen=20)

stratbr<-function(X,stratas=3,nsample=30,tampop=100,totgen=1500,pelite=0.2,
                  pmutant=0.3,rc=0.6,minimumsample=2,cores=2)
{
read_procedure<-function(X,Xmin)
  {X=sort(X)
    X=X[which(X>Xmin)]
    if (log10(max(X))>10000) {scale_div<<-max(X)} else {scale_div=1}
    X=X/scale_div
    B=unique(X);
    return(list(N=length(X),X=X,B=B,scale_mult=scale_div));
  }

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


  pop_gen<-function(psize,L,faixa)
  {
    popa<-cbind(matrix(runif(psize*(L-1),0,1),nrow=psize,ncol=L-1,byrow=TRUE),sample(faixa,psize,TRUE))
    return(popa)
  }


  decoder<-function(pop,stratas,Nb,B,X,qtd)
  {
    calcula_Vh_Nh<-function(x,X,H)
    {Vh<-rep(0,H)
    Nh<-rep(0,H)
    for(i in 1:H)
    {Xh<-which(X>x[i] & X<=x[i+1])
    Vh[i]<-var(X[Xh])
    if (is.na(Vh[i])==TRUE) {Vh[i]<-Inf}
    Nh[i]<-length(Xh)
    }
    return(cbind(Nh,Vh))
    }
    ak<-B[1]
    decod<-t(apply(pop,1,function(x) sort(ak+x[1:(stratas-1)]*(qtd[x[stratas]]-ak))))
    decod<-cbind(-Inf,decod,Inf)
    Nh_Vh<-matrix(apply(decod,1,function(x) calcula_Vh_Nh(x,X,stratas)),ncol=2*stratas,byrow=TRUE)
    return(cbind(Nh_Vh,decod))
  }



  ##########Objective Function#########################################################################
  fitness<-function(g,nsample,LL,minimumsample,X)
  {
    db<-t(snowfall::sfApply(g,1,function(g) optimal_allocation(g,nsample,LL,minimumsample,X)))
    return(db)
  }

  #########################################################################################################


  crossover_uniforme<-function(gelite,gnonelite,rc,pelite,pnelite,psort,N)
  {
    uniform<-function(ge,gn,ab,rc,N)
    {
      pr<-runif(N,0,1)
      gp<-which(pr<rc)
      ngp<-setdiff((1:N),gp)
      gn<-rep(0,N)
      gn[gp]<-gelite[ab[1],gp]
      gn[ngp]<-gnonelite[ab[2],ngp]
      return(gn)

    }

    ab=cbind(sample(pelite,psort,replace=TRUE),sample(pnelite,psort,replace=TRUE))
    gnew=t(apply(ab,1,function(ab) uniform(gelite,gnonelite,ab,rc,N)))
    return(gnew)
  }



  ##################Principal##########################
  snowfall::sfInit(parallel=TRUE,cpus=cores)
  population<-read_procedure(X,0)
  #N = comprimento do cromossomo = tamanho da populacao
  N<-length(population$B)
  tempo<-proc.time()
  size_population<-tampop
  pelite<-round(size_population*pelite)
  pmutant<-round(size_population*pmutant)
  pquantil<-1/length(population$B)
  pquantil<-ifelse(length(population$B)>100,0.01,0.1)
  if (length(population$B)>100) {faixa<-c(90:101)} else {faixa<-c(5:10)}
  qtd<-quantile(population$B,probs=seq(0,1,pquantil))
  f<-pop_gen(psize=size_population,stratas,faixa)
  tx<-sum(population$X)
  g<-decoder(f,stratas,N,population$B,population$X,qtd)
  ft<-fitness(g[,1:(2*stratas)],nsample,stratas,minimumsample,tx)
  fbest<-Inf
  i<-0
  smelhoriag<-0
  while((i<=totgen) & (smelhoriag<500))
     {i<-i+1
      smelhoriag<-smelhoriag+1
      pq<-order(ft[,(stratas+1)])
      f<-f[pq,] #Ordena pela Fitness
      fmin<-ft[pq[1],(stratas+1)]
      ft<-ft[pq,]
      g<-g[pq,]
      if (fmin<fbest)
       {fbest<-fmin
        gbest<-g[1,]
        solution<-ft[1,]
        cat("Solution Generation ",i," = ",round(fbest,6),"\n")
        flush.console()
        smelhoriag<-0
        ibest<-i
       }
  felite<-f[(1:pelite),] #Elite
  fnonelite<-f[(pelite+1):size_population,] #N?o Elite
  fmutant<-pop_gen(psize=pmutant,stratas,faixa) #Mutantes
  fnovos<-crossover_uniforme(felite,fnonelite,rc,pelite,(size_population-pelite),(size_population-pelite-pmutant),stratas)
  f<-rbind(fmutant,fnovos)
  g<-rbind(g[1:pelite,],decoder(f,stratas,N,population$B,population$X,qtd))
  f<-rbind(felite,f)
  glk<-g[(pelite+1):size_population,]
  ftk<-fitness(glk[,1:(2*stratas)],nsample,stratas,minimumsample,tx)
  ft<-rbind(ft[1:pelite,],ftk)

  }
  tempo<-(proc.time()-tempo)[3]
  snowfall::sfStop()
  return(list(cvtot=round(100*fbest,5),nh=solution[1:stratas],Nh=gbest[1:stratas],Sh2=gbest[(stratas+1):(2*stratas)],bk=gbest[(2*stratas+2):(3*stratas)],cputime=tempo))

}
