


### Clustering ultrasonic waves propagation time a hierarchical polynomial semi-parametric approach
This work is the result of a partnership between the Federal University of Ceará and the Federal University of São Calos,  Brazil, which experiment was carried out at the Laboratory of the Reabilitation and Durability of Constructions (LAREB),  <https://lareb.ufc.br/>.



Autors: Daiane Aparecida Zuanetti, Rosineide da Paz,  Talisson Rodrigues  and Esequiel Mesquita.






```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


## __Analysis of contamined simulated data set.__

In this work, a data-driven hierarchical regression model is applied to simulated data set. 



#### R packages required.

```{r echo=TRUE, paged.print=TRUE}
if(!require(mvtnorm)) install.packages("mvtnorm");require(mvtnorm) 
if(!require(MCMCpack)) install.packages("MCMCpack");require(MCMCpack) 
```


#### Definition of the simulated global population

```{r}
m=300    #number of individuals
Tp=5     #number of replications for each individual
n=m*Tp   #total number of observations
K=4      #number of groups with different effects
sem=100  ## seed
pcov<-2  # Total of fixed effects
q<-pran<-3  # Total of random effects
indiv=sort(rep(1:m,Tp)) # Information about individuals
repl=table(indiv)       # Information about replication in the individuals
#

```



#### Function for simulating indicator variables from a discrete distribution

```{r}
rDiscreta<-function(p){
 u<-runif(1)
 P<-cumsum(p)
 val<-sum(P<u)+1
 val}
#
pgrupos<-c(0.20,0.40,0.20,0.20)   ## Weight or probability for each group
set.seed(sem)
Sverd<-numeric(m)
for (i in 1:m){
	Sverd[i]<-rDiscreta(pgrupos)   ## Indicator group variable
}


```



#### Parameters for the fixed-effect

```{r}
set.seed(sem)
bfixo_verd =matrix(rnorm(pcov,0,1),ncol=1)
beta=c("beta1","beta2")
row.names(bfixo_verd)=beta


```



#### Parameters for the random-effect.


```{r}
a1=c(-4,-2,1,6)
a2=c(2,1.5,5,2)
set.seed(sem) 
sigm=diag(a2[1],q,q)
b1=rmvnorm(1,rep(a1[1],pran),sigm)
sigm=diag(a2[2],q,q)
b2=rmvnorm(1,rep(a1[2],pran),sigm)
sigm=diag(a2[3],q,q)
b3=rmvnorm(1,rep(a1[3],pran),sigm)
sigm=diag(a2[4],q,q)
b4=rmvnorm(1,rep(a1[4],pran),sigm)
#
Buv=t(rbind(b1,b2,b3,b4)) ## unique random-effect matrix
B<-NULL
for (i in 1:m){
	b<-Buv[,Sverd[i]]
	B=rbind(B,b)
}
```






#### Explained  variables for fixed-effect for each individual




```{r}
set.seed(sem)
X<-matrix(0,nrow=m,ncol=pcov)
X[,1]<-runif(m,0,1)
X[,2]<-rnorm(m)
X1f<-NULL
X2f<-NULL
for (i in 1:m){
  X1f<-c(X1f,rep(X[i,1],repl[i]))
  X2f<-c(X2f,rep(X[i,2],repl[i]))}
X<-cbind(X1f,X2f)
```



#### Simulated non-observed random variable for each group

## Random variable represent the effect of the group para Simulated non-observed random variables for each group, which represent the effects of the group.

```{r}
Vt<-1:Tp
aux=matrix(round(poly(Vt, degree=3),3),ncol=3)
aux<-cbind(1,aux)
Z<-NULL
for (i in 1:m) Z<-rbind(Z,aux)
Z=Z[,-1]
head(X); head(Z)
```





#### Generated response variable


```{r}
Y=NULL
set.seed(sem)
err=matrix(rnorm(n,0,1),m,Tp)
for(i in 1:m){
 	Y<-rbind(Y,(as.matrix(X[which(indiv==i),],nrow=sum(indiv==i))%*%bfixo_verd+
 	              Z[which(indiv==i),]%*%matrix(Buv[,Sverd[i]],ncol=1)+err[i,]) )}
#
nout=0  # if nout > 0 outliers are included in the sample
if(nout==0){sout=NULL}else{sout=1:nout
out=numeric()
set.seed(sem)
for(r in 1:nout){	out=rbind(out,(as.matrix(X[which(indiv==i),],nrow=sum(indiv==i))%*%bfixo_verd+
	              log(c(1,2,3,4,5))*(c(1,3,6,3,-8))  + rnorm(5,0,5)) )}
 Y[1:25]=out
}
#
ini=nout+1
plot(NULL, xlim=c(1,Tp), ylim=c(-15,15), ylab="Response", xlab="Index")
for (i in sout) lines(1:Tp,Y[indiv==i],col=(max(Sverd)+2),lty=(max(Sverd)+2))
for (i in ini:m) lines(1:Tp,Y[indiv==i],col=Sverd[i],lty=Sverd[i])
legend(1, 13.5, box.lty=0,legend=c("K = 1", "K = 2", "K = 3", "K = 4", "Outlier"),
       col=c(unique(Sverd),(max(Sverd)+2)),lty=c(unique(Sverd),(max(Sverd)+2)), cex=0.8,horiz=TRUE)
#
```




#### Function for sampling from the posterior of parameter $\sigma^2$

```{r}
poster_sigma2<-function(neta.a,neta.b,residuos){
 alpha<-(length(residuos)/2)+neta.a
 beta<-(sum(residuos^2)/2)+neta.b
 sigma2<-1/(rgamma(1,alpha,beta))
 return(sigma2)}
```


#### Function for sampling fro posterior distribution of parameter $\alpha$.

```{r}
gera_eta_alpha<-function(alpha,a_alph_prior,b_alph_prior,K,m){
	eta_aux<-rbeta(1,alpha+1,m)
	aux_prob<-(a_alph_prior+K-1)/(m*(b_alph_prior-log(eta_aux)))
	prob_alpha<-aux_prob/(1+aux_prob)
	unif<-runif(1)
	if (unif<=prob_alpha) alpha<-rgamma(1,a_alph_prior+K,b_alph_prior-log(eta_aux)) else alpha<-rgamma(1,a_alph_prior+K-1,b_alph_prior-log(eta_aux))
	return(alpha)}
```



##### Hyperparameters of the model

```{r}
lambda1<-lambda2<-0.01
qcov<-ncol(Z)
ni<-4 # degree of freedom of inverse-wishart
sigma<-diag(100,ncol=ncol(Z),nrow=ncol(Z)) #  variable and explained variable of inverse-wishart
for (i in 1:ncol(Z)){
	for (j in 1:ncol(Z)) if (sigma[i,j] != 100) sigma[i,j]<-0.5}

a_alph_prior<-b_alph_prior<-1
```


#### Useful informations for MCMC iterations.

```{r}
XXinv<-solve(t(X)%*%X)
Xtrans<-t(X)
XXinvXtrans<-XXinv%*%Xtrans
```


#### Initial values for  the MCMC chain



```{r}
library(mvtnorm)
library(MCMCpack)
#
set.seed(100)
S_vig<-rep(1,m)
S_obs_vig<-NULL
for (i in 1:m) S_obs_vig<-c(S_obs_vig,rep(S_vig[i],repl[i]))
K<-length(table(S_vig))
sigma2_vig<-round(poster_sigma2(lambda1,lambda2,c(Y)),3)
beta_vig<-matrix(0,ncol=1,nrow=qcov)
Dmat<-riwish(ni,sigma)
b_vig<-matrix(0,ncol=K,nrow=qcov)
for (k in 1:K){
	Bk<-Z[which(S_obs_vig==k),]
	aux<-solve((t(Bk)%*%Bk)+(solve(Dmat)*sigma2_vig))
	media<-aux%*%((t(Bk)%*%Y[which(S_obs_vig==k),])+(sigma2_vig*solve(Dmat)%*%beta_vig))
	varcov<-aux*sigma2_vig
	b_vig[,k]<-t(rmvnorm(1, mean = c(media), sigma = varcov, method=c("chol"), pre0.9_9994 = TRUE))}
ran_pred<-NULL
for (i in 1:n) ran_pred[i]<-Z[i,]%*%matrix(b_vig[,S_obs_vig[i]],ncol=1)
residuos<-matrix(Y-ran_pred,ncol=1)
gama_vig<-rmvnorm(1,mean = c(XXinvXtrans%*%residuos), sigma=(XXinv*sigma2_vig), method=c("chol"), pre0.9_9994 = TRUE)
alpha_vig<-gera_eta_alpha(1,a_alph_prior,b_alph_prior,K,m)
#
```



#### Definition for MCMC convergency.

```{r}
burnin<-500         ## burnig
amostrasfin<-1500    ## sample used for inference
saltos<-3            ## jump definition    
```


#### MCMC function for MCMC runing



```{r echo=TRUE, message=FALSE, warning=FALSE}
est.param <- list(sigma2=NULL, 
                  bs=list(),
                  Sj=NULL, 
                  Gamma=NULL, 
                  K=NULL, 
                  alpha=NULL, 
                  Betas=NULL, 
                  matrix_D=NULL,
                  log_vero=NULL)		
AmostrasTotal<-burnin+amostrasfin*saltos
K=1

set.seed(sem)

for (int in 1:AmostrasTotal){
#for (int in 1:8){
	cat('\n',c(int,K))
	#
	######## update S
	#
	for (i in 1:m){
		nk<-numeric(K)
		prob<-numeric(K)
		for (k in 1:K){
			nk[k]<-sum(S_vig[-i]==k)
			prob[k]<-nk[k]*dmvnorm(Y[which(indiv==i),],mean=c(matrix(X[which(indiv==i),],ncol=ncol(X))%*%matrix(gama_vig,ncol=1)+Z[which(indiv==i),]%*%matrix(b_vig[,k],ncol=1)),sigma=diag(sigma2_vig,repl[i]))}
		prob<-c(prob,alpha_vig*dmvnorm(Y[which(indiv==i),],mean=c(matrix(X[which(indiv==i),],ncol=ncol(X))%*%matrix(gama_vig,ncol=1)+Z[which(indiv==i),]%*%beta_vig),sigma=(Z[which(indiv==i),]%*%Dmat%*%t(Z[which(indiv==i),])+diag(sigma2_vig,repl[i]))))
		S_old<-S_vig[i]
		S_vig[i]<-rDiscreta(prob/sum(prob))
		if (S_vig[i]!=S_old){
			S_obs_vig[which(indiv==i)]<-S_vig[i]
			b_vig_teste<-matrix(0,nrow=qcov,ncol=max(K,max(S_vig)))
			b_vig_teste[,1:ncol(b_vig)]<-b_vig
			b_vig<-b_vig_teste
			for (k in 1:max(K,max(S_vig))){
				if (length(which(S_obs_vig==k))>0){
					Bk<-Z[which(S_obs_vig==k),]
					aux<-solve((t(Bk)%*%Bk)+(solve(Dmat)*sigma2_vig))
					media<-aux%*%((t(Bk)%*%(Y[which(S_obs_vig==k),]-matrix(X[which(S_obs_vig==k),],ncol=ncol(X))%*%matrix(gama_vig,ncol=1)))+(sigma2_vig*solve(Dmat)%*%beta_vig))
					varcov<-aux*sigma2_vig
					b_vig[,k]<-t(rmvnorm(1, mean = c(media), sigma = varcov, method=c("chol"), pre0.9_9994 = TRUE))}}}
		while (length(table(S_vig))<max(S_vig)){ # exclude empty clusters
			categr<-as.numeric(as.character(data.frame(table(S_vig))[,1]))
			categd<-seq(1:length(table(S_vig)))
			dif<-which(categr!=categd)
			S_vig[which(S_vig>dif[1])]<-S_vig[which(S_vig>dif[1])]-1
			b_vig<-matrix(c(b_vig[,-dif[1]]),nrow=qcov)
			K<-ncol(b_vig)
			S_obs_vig<-NULL
			for (ret in 1:m) S_obs_vig<-c(S_obs_vig,rep(S_vig[ret],repl[ret]))}	
		if (length(table(S_vig))<K){
			b_vig<-matrix(c(b_vig[,-K]),nrow=qcov)
			K<-ncol(b_vig)}
		K<-length(table(S_vig))}		
	#
    ###### ordena os grupos do maior para o menor
    #
    ord=as.numeric(names(sort(table(S_vig),decreasing = TRUE)))
  	Ss=S_vig
  	Ss_obs<-S_obs_vig
  	for(l in 1:length(ord)){
  		S_vig[Ss==ord[l]]=l
  		S_obs_vig[Ss_obs==ord[l]]=l}
	b_vig<-matrix(b_vig[,ord],ncol=ncol(b_vig))
    #
	###### update matrix D
	#
	Spost<-matrix(0,ncol=ncol(sigma),nrow=nrow(sigma))
	for (k in 1:K) Spost<-((matrix(b_vig[,k],ncol=1)-beta_vig)%*%t(b_vig[,k]-beta_vig))*sum(S_vig==k)+Spost
	Dmat<-riwish(ni+m,sigma+Spost)
	#
	###### update random effects
	#
	for (k in 1:K){
		Bk<-Z[which(S_obs_vig==k),]
		aux<-solve((t(Bk)%*%Bk)+(solve(Dmat)*sigma2_vig))
		media<-aux%*%((t(Bk)%*%(Y[which(S_obs_vig==k),]-matrix(X[which(S_obs_vig==k),],ncol=ncol(X))%*%matrix(gama_vig,ncol=1)))+(sigma2_vig*solve(Dmat)%*%beta_vig))
		varcov<-aux*sigma2_vig
		b_vig[,k]<-t(rmvnorm(1, mean = c(media), sigma = varcov, method=c("chol"), pre0.9_9994 = TRUE))}		
	#
	###### update fixed effects, variance and alpha
	#
	ran_pred<-NULL
	for (i in 1:n) ran_pred[i]<-Z[i,]%*%matrix(b_vig[,S_obs_vig[i]],ncol=1)
	media_beta<-matrix(0,ncol=1,nrow=nrow(b_vig))
	for (k in 1:K) media_beta<-media_beta+(sum(S_vig==k)*b_vig[,k])
	beta_vig<-t(rmvnorm(1, mean = c(media_beta/m), sigma = Dmat/m, pre0.9_9994 = TRUE))
	residuos<-matrix(Y-ran_pred,ncol=1)
	gama_vig<-rmvnorm(1,mean=c(XXinvXtrans%*%residuos), sigma=(XXinv*sigma2_vig), pre0.9_9994 = TRUE)
	residuos<-matrix(residuos-(X%*%matrix(gama_vig,ncol=1)),ncol=1)
	sigma2_vig<-round(poster_sigma2(lambda1,lambda2,c(residuos)),7)
	alpha_vig<-gera_eta_alpha(1,a_alph_prior,b_alph_prior,K,m)
	#
	#### recording results
	#	
	if (int>burnin & int%%saltos==0){
		log_vero<-0
		for (i in 1:m) log_vero<-log_vero+dmvnorm(Y[which(indiv==i),],mean=c(matrix(X[which(indiv==i),],ncol=ncol(X))%*%matrix(gama_vig,ncol=1)+Z[which(indiv==i),]%*%matrix(b_vig[,S_vig[i]],ncol=1)),sigma=diag(sigma2_vig,repl[i]),log=TRUE)
#		
		
b.list<-list(c(b_vig))
b.sample<-est.param[["bs"]]
est.param[["bs"]] <- c(b.sample,b.list)
est.param[["Sj"]]=rbind(est.param[["Sj"]],  S_vig)
est.param[["Gamma"]]=rbind(est.param[["Gamma"]],  gama_vig)
est.param[["K"]]=rbind(est.param[["K"]],  K)
est.param[["sigma2"]]=rbind(est.param[["sigma2"]],  sigma2_vig)
est.param[["alpha"]]=rbind(est.param[["alpha"]],  alpha_vig)
est.param[["Betas"]]=rbind(est.param[["Betas"]],  beta_vig)
est.param[["matrix_D"]]=rbind(est.param[["matrix_D"]], round(Dmat,3))
est.param[["log_vero"]]=rbind(est.param[["log_vero"]],  log_vero) 
}}
```








```{r}
iter<-amostrasfin
qcov<-ncol(Z)
pfixo<-ncol(X)
#
K<-est.param$K
Sj<-est.param$Sj
variancia<-est.param$sigma2
alfa<-est.param$alpha
betas<-est.param$Betas
baleat<-est.param$bs
mat_D<-est.param$matrix_D
gama<-est.param$Gamma
log_vero<-est.param$log_vero
```


#### Checking the convergence of the chain using the log-likelihood

```{r}
if(!require(coda)) install.packages("coda");require(coda) 
```




```{r}
log_v<-mcmc(log_vero)
effectiveSize(log_v)
geweke.diag(log_v)
plot(log_v,type='l')
```



#### The $alpha$'s estimates

```{r}
alfa<-mcmc(alfa)
effectiveSize(alfa)
geweke.diag(alfa)
plot(alfa,type='l')
mean(alfa)
quantile(alfa,c(0.025,0.50,0.975))
```


#### The $sigma^2$'s estimates


```{r}
variancia<-mcmc(variancia)
effectiveSize(variancia)
geweke.diag(variancia)
plot(variancia,type='l')
mean(variancia)
quantile(variancia,c(0.025,0.50,0.975))
```








#### Remaining parameters

```{r}
S<-matrix(Sj,ncol=m,nrow=iter,byrow=TRUE)
b_fixo<-matrix(0,nrow=iter,ncol=qcov)
b_aleat<-matrix(0,nrow=iter,ncol=(max(K)*qcov))
matrizD<-matrix(0,nrow=iter,ncol=qcov*qcov)
gama_v<-matrix(0,nrow=iter,ncol=pfixo)
obs<-1
obs2<-1
obs3<-1
for (i in 1:iter){
  if (K[i]>0){
    b_fixo[i,]<-betas[obs:(obs+qcov-1)]
    b_aleat[i,1:(K[i]*qcov)]<-baleat[[i]]
    matrizD[i,]<-c(t(mat_D))[(qcov*qcov*(i-1)+1):(qcov*qcov*(i))]
    gama_v[i,]<-gama[i,]}
  obs<-obs+qcov
  obs2<-obs2+(K[i]*qcov)
  obs3<-obs3+pfixo}
```









```{r}

#
####### fixed-effects' estimates
#
par(mfrow=c(2,1))
plot(gama_v[,1],type='l')
plot(gama_v[,2],type='l')
apply(gama_v,2,mean)
quantile(gama_v[,1],c(0.025,0.50,0.975))
quantile(gama_v[,2],c(0.025,0.50,0.975))
#
####### mean vector of G0's estimates
#
par(mfrow=c(3,1))
plot(b_fixo[,1],type='l')
plot(b_fixo[,2],type='l')
plot(b_fixo[,3],type='l')
apply(b_fixo,2,mean)
quantile(b_fixo[,1],c(0.025,0.50,0.975))
quantile(b_fixo[,2],c(0.025,0.50,0.975))
quantile(b_fixo[,3],c(0.025,0.50,0.975))
#
####### covariance matrix of G0's estimates
#
apply(matrizD,2,mean)
quantile(matrizD[,1],c(0.025,0.50,0.975))
quantile(matrizD[,2],c(0.025,0.50,0.975))
quantile(matrizD[,3],c(0.025,0.50,0.975))
quantile(matrizD[,4],c(0.025,0.50,0.975))
quantile(matrizD[,5],c(0.025,0.50,0.975))
quantile(matrizD[,6],c(0.025,0.50,0.975))
quantile(matrizD[,7],c(0.025,0.50,0.975))
quantile(matrizD[,8],c(0.025,0.50,0.975))
quantile(matrizD[,9],c(0.025,0.50,0.975))

```











#### Convergence analysis.

```{r}
log_vero <- est.param$log_vero 
log_v<-mcmc(log_vero)
effectiveSize(log_v)
geweke.diag(log_v)
```




#### Point estimation for parameter of the model, and indicator variable.

```{r}
Sigest <- est.param$sigma2 

mean(Sigest)
plot(Sigest,type="l")

Alpest <- est.param$alpha  
mean(Alpest)
plot(Alpest,type="l")

gamest <- est.param$Gamma  
colMeans(gamest)
plot(gamest[,1],type="l")

mean(Alpest)
plot(Alpest,type="l")

K <- est.param$K  
baleat <- est.param$bs 
Sj <- est.param$Sj 
betas <- est.param$Betas  
mat_D <- est.param$matrix_D   
iter<-length(K)

S<-Sj
b_fixo<-matrix(0,nrow=iter,ncol=pcov)
b_aleat<-matrix(0,nrow=iter,ncol=(max(K)*q))
matrizD<-matrix(0,nrow=iter,ncol=q*q)
obs<-1
obs2<-1
for (i in 1:iter){
  if (K[i]>0){
    b_fixo[i,]<-betas[obs:(obs+pcov-1)]
    b_aleat[i,1:(K[i]*q)]<-baleat[[i]]
    matrizD[i,]<-mat_D[(q*q*(i-1)+1):(q*q*(i))]}
  obs<-obs+q
  obs2<-obs2+(K[i]*q)}
#
Sj.j<-S # matrix of sample of indicator variable
prob.eq<-matrix(0,nrow=ncol(Sj.j),ncol=ncol(Sj.j))
for (i in 1:ncol(Sj.j)){
	for (j in 1:ncol(Sj.j)){
		prob.eq[i,j]<-round(sum(Sj.j[,i]==Sj.j[,j])/length(Sj.j[,i]),5)*100}}
#
thresh<-0.50*100 # set final groups
clust_f<-c(1,rep(0,(ncol(Sj.j)-1)))
for (i in 2:ncol(Sj.j)){
#for (i in 310:514){
	if (max(prob.eq[i,1:(i-1)])>thresh) clust_f[i]<-clust_f[which(prob.eq[i,1:(i-1)]==max(prob.eq[i,1:(i-1)]))[1]] else clust_f[i]<-max(clust_f[1:(i-1)]+1)}
table(clust_f) ###### final clusters
table(clust_f,Sverd)
#
#
###### For estimate random-effect of the groups
#
baleat_f<-matrix(0,nrow=iter,ncol=length(table(clust_f))*q)
for (gr in 1:length(table(clust_f))){
	ind_gr<-which(clust_f==gr)
	for (it in 1:iter){
		grupos<-S[it,ind_gr]
		for (cov in 1:q){
			baleat_f[it,((gr-1)*q+cov)]<-mean(b_aleat[it,((grupos-1)*q+cov)])
}}}
#
```

#### Interval estimates for the random-effect variables.


```{r}
#Elaboration of summary
R=Bestk<-baleat_f
Est=numeric()
for(i in 1:ncol(R)){
  vetor=pp<-matrix(R[,i])
  vetor <- mcmc(vetor)
  hpd=HPDinterval(vetor,prob=0.95)
  int=hpd[1:2]
  est=c(mean(R[,i]),int)
  Est=rbind(Est,est)
}



```




#### Relative frequency of the classification individuals.


```{r}
groups=as.numeric(names(sort(table(K), decreasing = FALSE)))

ran=max(unique(Sverd))
Tab=matrix(data=0, nrow=ran,ncol=length(groups))
for(j in 1:ran){  ### j representa os quadrntes de 1 até 24
idd=which(Sverd==j)  
for(i in idd){
Tab[j,Se[i]]=Tab[j,Se[i]]+1
}
}

for(i in 1:nrow(Tab)){Tab[i,]=Tab[i,]/sum(Tab[i,])}

Tab
```




#### Defining final clusters

```{r}


#
####### Defining final clusters
#
Sj.j<-S #co-clustering posterior probabilities
prob.eq<-matrix(0,nrow=ncol(Sj.j),ncol=ncol(Sj.j))
for (i in 1:ncol(Sj.j)){
	for (j in 1:ncol(Sj.j)){
		prob.eq[i,j]<-round(sum(Sj.j[,i]==Sj.j[,j])/length(Sj.j[,i]),5)*100}}
#
thresh<-0.50*100 # thresholding co-clustering probability to cluster individuals
clust_f<-c(1,rep(0,(ncol(Sj.j)-1)))
for (i in 2:ncol(Sj.j)){
	if (max(prob.eq[i,1:(i-1)])>thresh) clust_f[i]<-clust_f[which(prob.eq[i,1:(i-1)]==max(prob.eq[i,1:(i-1)]))[1]] else clust_f[i]<-max(clust_f[1:(i-1)]+1)}
table(clust_f,Sverd)
#
thesing<-0.3 # merging atypical individuals into groups where at least 30% of MCMC iterations they appear
singl<-which(clust_f %in% which(table(clust_f)==1))
prob.eq.sin<-matrix(prob.eq[singl,],nrow=length(singl))
for (i in 1:nrow(prob.eq.sin)){
	prob.eq.sin[i,singl[i]]<-0
	if (max(prob.eq.sin[i,])>thesing) clust_f[singl[i]]<-clust_f[which(prob.eq.sin[i,]==max(prob.eq.sin[i,]))[1]]}
while (length(table(clust_f))<max(clust_f)){ # exclude empty clusters
	categr<-as.numeric(as.character(data.frame(table(clust_f))[,1]))
	categd<-seq(1:length(table(clust_f)))
	dif<-which(categr!=categd)
	clust_f[which(clust_f>dif[1])]<-clust_f[which(clust_f>dif[1])]-1}
#
table(clust_f,Sverd) # final clusters
#
####### Random-effects' estimates for each group
#
baleat_f<-matrix(0,nrow=iter,ncol=length(table(clust_f))*qcov)
for (gr in 1:length(table(clust_f))){
	ind_gr<-which(clust_f==gr)
	for (it in 1:iter){
		grupos<-S[it,ind_gr]
		for (cov in 1:qcov){
			baleat_f[it,((gr-1)*qcov+cov)]<-mean(b_aleat[it,((grupos-1)*qcov+cov)])
}}}
B_aleat<-matrix(apply(baleat_f,2,mean),
                ncol=length(table(clust_f)),nrow=qcov)
B_aleat
#

```




```{r}

#Elaboration of summary
R=Bestk<-baleat_f
Est=numeric()
for(i in 1:ncol(R)){
  vetor=pp<-matrix(R[,i])
  vetor <- mcmc(vetor)
  hpd=HPDinterval(vetor,prob=0.95)
  int=hpd[1:2]
  est=c(mean(R[,i]),int)
  Est=rbind(Est,est)
}



v=round(Est,2)
paste(v[1,1],"(",v[1,2],",",v[1,3],")","&",
      v[2,1],"(",v[2,2],",",v[2,3],")","&",
      v[3,1],"(",v[3,2],",",v[3,3],")","&",
      sep=" ")
```








#### Graphical representation of final estimated individuals classification.


```{r}

# label organization
clust_f2<-NULL

for (i in 1:length(clust_f)){

	if (clust_f[i]==1) clust_f2[i]<-2

	if (clust_f[i]==2) clust_f2[i]<-1

	if (clust_f[i]==3) clust_f2[i]<-4

	if (clust_f[i]==4) clust_f2[i]<-3		

}


#

baleat_f2<-matrix(0,nrow=nrow(baleat_f),ncol=ncol(baleat_f))

baleat_f2[,1]<-baleat_f[,2]

baleat_f2[,2]<-baleat_f[,1]

baleat_f2[,3]<-baleat_f[,4]

baleat_f2[,4]<-baleat_f[,3]




```






```{r}

Se=clust_f2
#jpeg("output.jpg", width = 1000, height = 10000)
#pdf("simulatedplot.pdf", width = 8, height = 5)
par(mar=c(4,5,0,0))
plot(NULL, xlim=c(1,Tp), ylim=c(-15,16), ylab="Response", xlab="Index",cex.axis=2,cex.lab=2)
for (i in 1:m) lines(1:Tp,Y[indiv==i],col=Se[i],lty=Se[i],lwd=unique(Se))
legend(1, 16.5, box.lty=0,legend=unique(Se),
       col=unique(Se),lty=unique(Se), cex=2,horiz=TRUE,lwd=unique(Se))
#dev.off()

#jpeg("output.jpg", width = 1000, height = 10000)
#pdf("realplot.pdf", width = 8, height = 5)
par(mar=c(4,5,0,0))
plot(NULL, xlim=c(1,Tp), ylim=c(-15,16), ylab="Response", xlab="Index",cex.axis=2,cex.lab=2)
for (i in 1:m) lines(1:Tp,Y[indiv==i],col=Sverd[i],lty=Se[i],lwd=unique(Sverd))
legend(1, 16.5, box.lty=0,legend=unique(Sverd),
       col=unique(Sverd),lty=unique(Sverd), cex=2,horiz=TRUE,lwd=unique(Sverd))
#dev.off()
```
