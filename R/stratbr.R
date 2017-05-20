#'Optimization Algorithm to solve stratification problem
#'
#'@param X             Stratification variable.
#'@param H             Number of strata.
#'@param n             Sample size.
#'@param nmin          Minimum sample size (smallest possible sample size in any stratum).
#'@param takeall       Take-all stratum (takeall=TRUE) =>  nH=NH.
#'@param tampop        Number of chromosomes BRKGA.
#'@param totgen        Number of generations BRKGA.
#'@param pelite        Percentage elite solutions BRKGA.
#'@param pmutant       Percentage mutant solutions BRKGA.
#'@param rc            Crossover probability  BRKGA.
#'@param cores         Numerical amount of CPUs requested for the cluster.
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
#'The optimization algorithm is applied to solve the one dimensional case,
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
#'solution1<-stratbr(REV84,H=3,n=50,nmin=10,totgen=2,cores=4)
#'data(USbanks)
#'solution2<-stratbr(USbanks,H=3,n=50,totgen=2,cores=4,takeall=TRUE)

stratbr<-function(X,H=3,n=30,nmin=2,takeall=FALSE,tampop=100,totgen=1500,pelite=0.2,
                  pmutant=0.3,rc=0.6,cores=2)
{
read_procedure<-function(X,Xmin)
  {X=sort(X)
    X=X[which(X>Xmin)]
    if (log10(max(X))>10000) {scale_div<<-max(X)} else {scale_div=1}
    X=X/scale_div
    B=unique(X);
    return(list(N=length(X),X=X,B=B,scale_mult=scale_div));
  }


pop_gen<-function(psize,L,faixa)
  {
    popa<-cbind(matrix(runif(psize*(L-1),0,1),nrow=psize,ncol=L-1,byrow=TRUE),sample(faixa,psize,TRUE))
    return(popa)
  }


  decoder<-function(pop,H,Nb,B,X,qtd)
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
    decod<-t(apply(pop,1,function(x) sort(ak+x[1:(H-1)]*(qtd[x[H]]-ak))))
    if (H==2) {decod<-t(decod)}
    decod<-cbind(-Inf,decod,Inf)
    Nh_Vh<-matrix(apply(decod,1,function(x) calcula_Vh_Nh(x,X,H)),ncol=2*H,byrow=TRUE)
    return(cbind(Nh_Vh,decod))
  }



  ##########Objective Function#########################################################################
  fitness<-function(g,n,LL,nmin,X,takeall)
  { db<-t(snowfall::sfApply(g,1,function(g) BSSM_FD(g[1:LL],g[(LL+1):(2*LL)],n,LL,nmin,X,takeall=takeall)))
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
  f<-pop_gen(psize=size_population,H,faixa)
  tx<-sum(population$X)
  g<-decoder(f,H,N,population$B,population$X,qtd)
  ft<-fitness(g[,1:(2*H)],n,H,nmin,tx,takeall)
  fbest<-Inf
  i<-0
  smelhoriag<-0
  while((i<=totgen) & (smelhoriag<500))
     {i<-i+1
      smelhoriag<-smelhoriag+1
      pq<-order(ft[,(H+1)])
      f<-f[pq,] #Ordena pela Fitness
      fmin<-ft[pq[1],(H+1)]
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
  fmutant<-pop_gen(psize=pmutant,H,faixa) #Mutantes
  fnovos<-crossover_uniforme(felite,fnonelite,rc,pelite,(size_population-pelite),(size_population-pelite-pmutant),H)
  f<-rbind(fmutant,fnovos)
  g<-rbind(g[1:pelite,],decoder(f,H,N,population$B,population$X,qtd))
  f<-rbind(felite,f)
  glk<-g[(pelite+1):size_population,]
  ftk<-fitness(glk[,1:(2*H)],n,H,nmin,tx,takeall)
  ft<-rbind(ft[1:pelite,],ftk)

  }
  tempo<-(proc.time()-tempo)[3]
  snowfall::sfStop()
  return(list(cvtot=round(100*fbest,5),nh=solution[1:H],Nh=gbest[1:H],Sh2=gbest[(H+1):(2*H)],bk=gbest[(2*H+2):(3*H)],cputime=tempo))
}


#'Optimal Allocation - Minimum Coefficient of Variation
#'
#'@param Nh            Vector with number of population elements, or population size, in stratum h
#'@param Sh2x          Vector with population variance of the variable X in stratum h.
#'@param n             Sample size.
#'@param H             Number of strata.
#'@param nmin          Minimum sample size (smallest possible sample size in any stratum).
#'@param X             Population Total
#'@param takeall       Take-all stratum (takeall=TRUE) =>  nH=NH.
#'@return \item{solution}{Vector with  sample of size by stratum and coefficient of variation
#'                        for the estimator of total of the stratification variable considered.}
#'
#'@references Brito, J.A.M, Silva, P.L.N.,Semaan, G.S. and Maculan, N. (2015).
#'            Integer Programming Formulations Applied to Optimal Allocation in Stratified Sampling.
#'            Survey Methodology, 41: 427-442.
#'@author Jose Brito (jambrito@gmail.com), Pedro Silva, Gustavo Semaan and Nelson Maculan.
#'
#'@description Function that uses an integer programming formulation for allocation
#'       of the overall sample size n to the strata,  for the following purpose:
#'    Coefficient of Variation of the estimate of total for the survey variable is minimized.
#'@export
#'@examples
#'X<-round(100*runif(50))
#'Nh<-c(10,20,20)
#'Sh2x<-c(var(X[1:10]),var(X[11:30]),var(X[31:50]))
#'aloc1<-BSSM_FD(Nh,Sh2x,n=40,H=3,nmin=2,sum(X),takeall=TRUE)
#'Nh<-c(49,78,20,39,73,82,89)
#'X<-542350
#'Sh2x<-c(4436978,5581445,33454902,5763294,8689167,3716130,13938505)
#'aloc2<-BSSM_FD(Nh,Sh2x,n=100,H=7,nmin=2,X)

BSSM_FD<-function(Nh,Sh2x,n,H,nmin,X,takeall=FALSE)
{
  Nh=Nh
  Vh=Sh2x
  funcao_objetivo=as.vector(outer(1/c(1:n),(Vh*Nh^2)))
  inf=1
  sup=n
  for(i in 1:H) {funcao_objetivo[inf:sup]=funcao_objetivo[inf:sup]-(Vh[i]*Nh[i]);inf=sup+1;sup=sup+n}
  funcao_objetivo[which(funcao_objetivo<0)]=max(Vh)*max(Nh^2)*10
  restricoes=c(rep("==",H+1),rep("<=",H),rep(">=",H))

  b=c(rep(1,H),n,Nh,rep(nmin,H))
  A=matrix(0,nrow=(3*H)+1,ncol=n*H)
  for(i in 1:H) {A[i,(1+(n*(i-1))):(i*n)]=rep(1,n)}
  A[H+1,]=rep(c(1:n),H)
  inf=1
  sup=n
  for(i in (H+2):((2*H)+1)){A[i,inf:sup]=c(1:n);inf=sup+1;sup=sup+n}
  for(i in ((2*H)+2):((3*H)+1)) {A[i,]=A[i-H,]}
  if (takeall==TRUE)
  {restricoes[2*H+1]<-"=="}

  solucao=Rglpk_solve_LP(funcao_objetivo,A,restricoes,b,types=rep("B",length(funcao_objetivo)),max=FALSE)
  status=solucao$status
  var_tot=solucao$optimum
  solucao=solucao$solution
  x=matrix(rep(0,H),nrow=1,byrow=TRUE)
  x[1,1]=sum(solucao[1:n]*c(1:n))
  if (H>1)
  {inf=n+1
  sup=2*n
  for(i in 2:H)
  {x[1,i]=sum(solucao[inf:sup]*c(1:n))
  inf=sup+1
  sup=sup+n
  }
  }

  if (status>0) {var_tot=Inf}
  solution=c(x,sqrt(var_tot)/X)
  return(solution)
}
