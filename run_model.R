

run_model <- function(data_list, MCMC=TRUE){
  ##
  if(isTRUE(MCMC)){
    stan_code <- "
      data {
        //
          int <lower = 0> N; // Defining the number of data points
          array[N] int<lower=0, upper=1> y; // A variable that describes [1] or [0] 
          int <lower = 0> K;   // number of predictors
          int <lower = 0> Q;   // number of cluster level predictors
          int <lower = 0> J;   // number of clusters
          int <lower = 0> M;   // number of 2nd stage areal mid-points
          matrix[N, K] x;   // K predictor matrix
          matrix[J, Q] z;   // Q predictors at J cluster/spatial scale 
          matrix[J, Q] zstar;   // Q predictors at J cluster/spatial scale 
          array[J] int<lower=1, upper=2> ur; // urban-rural ID; 1=urban, 2=rural
          matrix[J, M] D_JxM; // distance matrix
          matrix[J, M] Dstar; // distance matrix
          matrix[N, J] IMAT; // matrix
          //
          matrix[J, J] MJJ; // matrix
          matrix[M, Q] Sigma_diag; // matrix
          real phi; // decay 
          //
      }
    parameters {
      vector[K] beta;   // coefficients for predictors (including intercept, if applicable)
      vector[Q] zeta;   // coefficients for spatial predictors (including intercept, if applicable)
      real <lower = 0> sigma_nu_hat; // for zstar - true cluster level values
      real <lower = 0> sigma_u_hat; //  initiate - urban
      real <lower = 0> sigma_r_hat; //  initiate - rural
    }
    transformed parameters {
      vector[K] OR_beta;   // OR of coefficients for predictors (including intercept, if applicable)
      vector[Q] OR_zeta;   // OR of coefficients for spatial predictors (including intercept, if applicable)
      OR_beta = exp(beta);
      OR_zeta = exp(zeta);
    }
    model {
      //
      sigma_u_hat ~ inv_gamma(20, 0.05); // prior
      sigma_r_hat ~ inv_gamma(10, 0.04); // prior
      sigma_nu_hat ~ inv_gamma(2, 1); // prior
      //
      matrix[J, M] Dstar_phi; // matrix
      //
        for(m in 1:M){
          for(j in 1:J){
            if(ur[j] == 1){ // rural
              Dstar[j, m] ~ normal(D_JxM[j, m], sigma_r_hat);
            }
            else{ // urban
              Dstar[j, m] ~ normal(D_JxM[j, m], sigma_u_hat);
            }
            Dstar_phi[j, m] = (1-(Dstar[j, m]/phi)^2)^2;
          }
        }
      //
      matrix[J, M] Psi; 
      matrix[J, Q] Psi_eta; 
      matrix[J, Q] zstar_hat; 
      Psi = MJJ * Dstar_phi; // (JxJ) x (JxM) = JxM
      Psi_eta = Psi * Sigma_diag; // (JxM) x (MxQ) = JxQ
      //
        for(q in 1:Q){
          for(j in 1:J){
            zstar_hat[j, q] ~ normal(z[j, q] + Psi_eta[j, q], sigma_nu_hat); // prior for zstar 
            z[j, q] ~ normal(zstar_hat[j, q] - Psi_eta[j, q], sigma_nu_hat); 
          }
        }
      // 
      matrix[N, Q] zzstar; 
      zzstar = IMAT * zstar_hat; // (NxJ) x (JxQ) = NxQ // IMAT => 0 and 1 design matrix
      //
        for(i in 1:N){
          y[i] ~ bernoulli_logit(x[i,] * beta + zzstar[i,] * zeta); // likelihood
        }
      //
        for(k in 1:K){
          beta[k] ~ normal(0, 2);
        }
      for(q in 1:Q){
        zeta[q] ~ normal(0, 2);
      }
      //
    }
    generated quantities {
      vector[N] y_pred_prob;
      int<lower=0, upper=1> y_pred[N];
      vector[N] log_lik;
      vector[J] z_effect;
      vector[N] zzstar_hat;
      z_effect = zstar * zeta;          // Jxq * q => J
      zzstar_hat = IMAT * z_effect;     // NxJ * J => N
      for (i in 1:N) {
        real eta = dot_product(x[i], beta) + zzstar_hat[i];
        y_pred_prob[i] = inv_logit(eta);
        y_pred[i] = bernoulli_rng(y_pred_prob[i]);
        log_lik[i] = bernoulli_logit_lpmf(y[i] | eta);
      }
    }
	"
	out <- stan(
        model_code = stan_code,
        data = data_list,      # named list of data
        chains = 4,            # number of Markov chains
        warmup = 2000,         # number of warmup iterations per chain
        iter = 5000,           # total number of iterations per chain
        cores = 4              # number of cores (could use one per chain)
        #refresh = 0           # no progress shown
        )
	out	
  }
  else{
    stop("Currently working on this extension for the next paper ...")  
  }
  ##
}
