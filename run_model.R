

run_model <- function(model, data_list, digits=4, draws = 100, MCMC=FALSE){
  ##
  if(isTRUE(MCMC)){
    stop("Will add it ...")
  }
  else{
    K <- data_list$K
    Q <- data_list$Q
    out <- optimizing(model, data = data_list, hessian = TRUE, draws = draws)
    opt_cov <- solve(-out$hessian[1:(K+Q),1:(K+Q)])
    post.beta.mn <- out$par[1:K]
    post.beta.sd <- sqrt(diag(opt_cov[1:K,1:K]))
    parsim <- cbind(post.beta.mn,post.beta.sd)
    myNorm_exp <- function(x, sim=10000){
      quantile(exp(rnorm(sim,x[1],x[2])),prob=c(0.5,0.025,0.975))
    }
    para < round(t(apply(parsim,1,myNorm_exp)),2)
    dimnames(para)[[1]] = dimnames(data_list$x)[[2]]
    post.zeta.mn <- out$par[(K+1):(K+Q)]
    post.zeta.sd <- sqrt(opt_cov[(K+1):(K+Q),(K+1):(K+Q)])
    para_z <- quantile(exp(rnorm(1000, post.zeta.mn, post.zeta.sd)), prob = c(0.5, 0.025, 0.975))
    para_summary <- rbind(para,para_z)
    dimnames(para_summary)[[1]][nrow(para_summary)] <- cluster_vars_idz[2]
    round(para_summary,digits=digits)
  }
  ##
}
