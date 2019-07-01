#' @title Runs the integrative CAR Horseshoe model
#'
#' @description Infers treatment effects, association with heterogeneous omic variables, pathway perturbation
#' among other parameters (e.g. time dependence). Regression coefficients (beta parameter) are initialized
#' using a univariate regression ignoring time and metabolite dependence.
#'
#' @param X the metabolomics time-course data with dimensions timepoints x observations x variables
#' @param Y the additional omic time-course data with dimensions timepoints x observations x variables
#' @param drug treatment effect (NA values not allowed in drug) with dimensions timepoints x observations
#' @param groups grouping vector (binary).
#' @param pathways pathway adjacency matrices as returned by iCARH.getPathwaysMat
#' @param tau global sparsity parameter \eqn{\tau} as in Jendoubi, T., & Ebbels, T. (2018)
#' @param NA_value NA values are incompatible with stan.
#' NAs will be replaced by NA_value and will be inferred (only for X and Y data).
#' @param init If \code{TRUE} use iCARH provided initialization function. Passed to Stan otherwise. Please see Stan manual 
#' on \code{init} possible values.
#' @param ... additional stan parameters
#'
#' @return stan object
#'
#' @examples data.sim = iCARH.simulate(4, 8, 10, 2, 2, path.probs=0.3, Zgroupeff=c(0,4),
#' beta.val=c(1,-1,0.5, -0.5))
#' XX = data.sim$XX
#' Y = data.sim$Y
#' Z = data.sim$Z
#' pathways = data.sim$pathways
#' \donttest{
#' rstan_options(auto_write = TRUE)
#' options(mc.cores = 2)
#' fit = iCARH.model(XX, Y, Z, pathways, control = list(adapt_delta = 0.99, max_treedepth=10),
#' iter = 2, chains = 2)}
#'
#' @export iCARH.model
#' @importFrom rstan stan
#' @importFrom utils packageVersion

iCARH.model = function(X, Y, drug, groups=NULL, pathways, tau=1.2, NA_value=-99999, init=T, ...){

  if ((packageVersion("rstan") <= "2.18.2") & (getRversion() >= "3.6.0")) {
    warning("rstan > 2.18.2 or R < 3.6.0 needed for this function. Program will exit promptly.", call. = FALSE)
    return(NULL)
  }

  if(is.null(groups)){
    if(is.logical(drug) | setequal(drug, c(0,1))) groups=!drug[1,]
    else warning("Need group information in logical or binary form.")
  }
  
  Xraw=X
  Yraw=Y
  x_i_mis = which(is.na(X), arr.ind = T);
  x_i_obs = which(!is.na(X), arr.ind = T);
  X[which(is.na(X))] =  NA_value;
  y_i_mis = which(is.na(Y), arr.ind = T);
  y_i_obs = which(!is.na(Y), arr.ind = T);
  Y[which(is.na(Y))] =  NA_value;
  adjmat = lapply(pathways, function(x) 1/(x+1))
  adjmat = lapply(adjmat, function(x) {diag(x)=0; 1/max(rowSums(x>0))*x})
  lambdas = lapply(adjmat,function(x) sort(eigen(x)$values)[c(1,nrow(x))]);
  cat.groups = as.numeric(as.character(factor(groups, levels=sort(unique(groups), decreasing=T), labels=1:length(unique(groups)) )))
  regression.stan = "
  data {
  int<lower=0> J; // metabolites
  int<lower=0> T; // timepoints
  int<lower=0> N; // observations
  int<lower=0> P; // pathways
  int<lower=0> K; // Y variables
  int<lower=1> G; // Groups, limited to 2 groups currently
  // Data and missing data
  int<lower=0> x_n_mis;
  int<lower=0> y_n_mis;
  matrix[N,J]  X[T];
  int<lower=1> x_i_mis[x_n_mis,3];
  int<lower=1> y_i_mis[y_n_mis,3];
  matrix[N,K]  Y[T];
  vector[N]  drug[T];
  real  NA_value;
  int<lower=1> groups[N];
  // Pathway data
  matrix<lower=0>[J,J] adjmat[P];
  vector[2] lambdas[P];
  // Overall shrinkage parameter
  real<lower=1> nu;
  }
  transformed data{
  matrix[J,J] I;
  I = diag_matrix(rep_vector(1,J));
  }
  parameters {
  // missing data
  real Xmis[x_n_mis];
  real Ymis[y_n_mis];
  // model parameters of interest
  real<lower=-1, upper=1> theta[J];
  vector[K] beta_std[J];
  real alpha[J];
  real inter[J]; // fixed intercept
  // variances
  real<lower=0> sigma;
  real<lower=0> sigma_gamma[J];
  real<lower=0> sigma_y;
  vector[N] err[J];
  // pathway coeff
  vector<lower=0, upper=1>[2] phi_std[P];
  // shrinkage prior parameters
  vector<lower=0>[J] nu1_g;
  vector<lower=0>[J] nu2_g;
  vector<lower=0>[K] nu1_l[J];
  vector<lower=0>[K] nu2_l[J];
  }
  transformed parameters{
  matrix[N,J] XX[T];
  matrix[N,J] Xm[T];
  matrix[N,K] YY[T];
  matrix[N,J] mu[T];
  vector[N] gamma[J]; // random intercept
  vector[K] beta[J];
  vector[G] phi[P];
  cov_matrix[J] Sigma[G];
  vector<lower=0>[K] nu12_l[J];
  vector<lower=0>[J] sigma_beta;
  // dealing with missing data
  for (i in 1:T){
  for(j in 1:J) {for(n in 1:N) if((X[i,n,j])!= NA_value) XX[i,n,j]=X[i,n,j];}
  }
  for(i in 1:T){
  for(k in 1:K) {for(n in 1:N) if((Y[i,n,k])!= NA_value) YY[i,n,k]=Y[i,n,k];}
  }
  for(l in 1:x_n_mis){
  XX[x_i_mis[l,1], x_i_mis[l,2], x_i_mis[l,3]] = Xmis[l];
  }
  for(l in 1:y_n_mis){
  YY[y_i_mis[l,1], y_i_mis[l,2], y_i_mis[l,3]] = Ymis[l];
  }
  // covariance for CAR level
  for(g in 1:G) {
  Sigma[g] = rep_matrix(0,J,J);
  for(p in 1:P){
  phi[p,g] = 1/(lambdas[p,1])+0.005 + (1/(lambdas[p,2]) - 1/(lambdas[p,1])-0.005) * phi_std[p,g];
  Sigma[g] = Sigma[g] + phi[p,g]*adjmat[p];
  }
  Sigma[g] = (I-Sigma[g]/P)/sigma;
  }
  // shrinkage prior, global and local levels
  sigma_beta = nu1_g .* sqrt(nu2_g);
  for(j in 1:J){
  gamma[j] = inter[j] + sqrt(sigma_gamma[j])*err[j];
  nu12_l[j] = nu1_l[j] .* sqrt(nu2_l[j]);
  beta[j] = beta_std[j].*nu12_l[j] * sigma_beta[j];
  }
  // mean for CAR level
  for (i in 1:(T)){ for(j in 1:J){ Xm[i,1:N,j] = gamma[j] + YY[i]*(beta[j]) +alpha[j]*drug[i];}}
  mu[1] = Xm[1];
  for (i in 2:(T)){ for(j in 1:J){ mu[i,1:N,j] = Xm[i,1:N,j] + theta[j]*(XX[i-1,1:N,j]-Xm[i-1,1:N,j]);}}
  }
  model {
  //priors on pathway significance
  for(g in 1:G) { for(p in 1:P) phi_std[p,g] ~ beta(0.5, 0.5);}
  //variances
  sigma ~ inv_gamma(N*T/4.0,N*T/4.0-1);
  nu1_g ~ normal(0.0, 1.0);
  nu2_g ~ inv_gamma(0.5, 0.5);
  sigma_y ~ inv_gamma(1,1);
  sigma_gamma ~ inv_gamma(1,0.1);
  for(j in 1:J){
  err[j] ~ normal(0,1);
  nu1_l[j] ~ normal(0.0, 1.0);
  nu2_l[j] ~ inv_gamma(0.5*nu, 0.5*nu);
  beta_std[j] ~ normal(0, 1);
  }
  // Inferring missing values in YY
  for(k in 1:K) YY[1,,k] ~ multi_normal(rep_vector(0,N),diag_matrix(rep_vector(1,N)));
  for(i in 2:T){
  for(k in 1:K) target += multi_normal_lpdf(YY[i,,k] | YY[i-1,,k], diag_matrix(rep_vector(sqrt(sigma_y),N)));
  }
  // Inference
  for (i in 1:(T)){
  for(n in 1:N){
  target += multi_normal_prec_lpdf(XX[i,n] | mu[i,n], Sigma[groups[n]]);
  }
  }
  }
  generated quantities {
  // checks
  real log_lik[T,N];
  matrix[J,N] mupred[T];
  for (i in 1:T){
  for(n in 1:N) log_lik[i,n] = multi_normal_prec_lpdf(XX[i,n] | mu[i,n], Sigma[groups[n]]);
  for(n in 1:N) mupred[i,,n] = multi_normal_rng((mu[i,n])', Sigma[groups[n]]);
  }
  }"
  if (init==T){
  initf = function(){
    #if(perturb)
    X[which(X== NA_value)] = NA
    Y[which(Y== NA_value)] = NA
    X = X + array(rnorm(prod(dim(X)),0,sd(X,na.rm=T)/10), dim=dim(X))
    coeff = array(dim=c(4,dim(X)[3], dim(Y)[3]))
    for(k in 1:dim(Y)[3]){
    for(i in 1:dim(X)[3]){
      dd = lm(X~.,data=data.frame(X=as.vector(X[2:nrow(X),,i]), as.vector(X[1:(nrow(X)-1),,i]),
                                    as.vector(Y[1:(nrow(X)-1),,k]), as.vector(drug[1:(nrow(X)-1),])))
      coeff[,i,k] = dd$coefficients
    }
    }
    return(list(inter=rowMeans(coeff[1,,]),theta=rowMeans(coeff[2,,]),beta=coeff[3,,],alpha=rowMeans(coeff[4,,])))
  }
  } else{
    initf=init
  }
dat = list(N=ncol(X), T=nrow(X), J=dim(X)[3], P=length(pathways), K=dim(Y)[3], G=max(cat.groups), X=X, Y=Y,
           drug=drug, groups=cat.groups, lambdas=lambdas, adjmat=adjmat, x_i_mis=x_i_mis, y_i_mis=y_i_mis,
           x_n_mis=nrow(x_i_mis), y_n_mis=nrow(y_i_mis), nu=tau,  NA_value= NA_value)
fit = rstan::stan(model_name="iCARH",model_code = regression.stan, data = dat, init=initf, pars=c("Xmis","Ymis"), include=F, ...)
fit = list(icarh=fit, X=Xraw, Y=Yraw, drug=drug, groups=groups)
class(fit) = "icarh"
return(fit)
  }
