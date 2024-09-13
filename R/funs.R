
#' Logit Transformation
#' This function computes the logit transformation of a given numeric input.
#' @export
logit <- function(x) {log(x/(1-x))}

# Function to do a perturbation for a given subgroup
adjust_responses <- function(n, r, p=0.5) {
  if (r == 0) { return(r + 1) }
  if (r == n) { return(r - 1) }
  if (runif(1) < p) { return(r + 1) } else { return(r - 1) }
}

# Function to calculate kullback leibler divergence
kld <- function(p, q, base = c("log", "log2", "log10"), margin = FALSE) {
  base <- match.arg(base)
  klds <- p * do.call(base, list(p / q))
  if (margin) {
    return(sum(klds, na.rm = TRUE))
  } else {
    klds
  }
}

# Function to calculate jesson shannon divergence
jsd <- function(p, q, base = c("log", "log2", "log10"), margin = FALSE) {
  base <- match.arg(base)
  m <- (p + q) / 2
  jsds <- 0.5 * (kld(p , m, base) + kld(q , m, base))
  if (margin) {
    return(sum(jsds, na.rm = TRUE))
  } else {
    jsds
  }
}

# Function to calculate jesson shannon divergence between two beta distribution
jsd.two.beta <- function(n, r1, r2, aprior=0.5, bprior=0.5){
  length <- 0.001
  interval <- seq(0, 1, by = length)
  x <- interval[-1]
  aposterior1 <- aprior + r1
  bposterior1 <- bprior + n - r1
  aposterior2 <- aprior + r2
  bposterior2 <- bprior + n - r2
  dist1 <- dbeta(x, aposterior1, bposterior1)*length
  dist2 <- dbeta(x, aposterior2, bposterior2)*length
  jsd(dist1, dist2, base = "log2", margin = T)
}

#' Function to filter extremely small groups
#'
#' This function filters small groups based on a perturbation analysis
#' and Jensen-Shannon Divergence (JSD) using the DBSCAN clustering algorithm.
#'
#' @param group A vector of subgroup labels.
#' @param n A vector representing the number of patients in each subgroup.
#' @param r A vector representing the number of responses in each subgroup.
#' @param adj_p Numeric. Adjustment parameter for adjusting responses.
#' @param aprior Numeric. Prior parameter 'a' for Bayesian estimation.
#' @param bprior Numeric. Prior parameter 'b' for Bayesian estimation.
#' @param eps Numeric. The epsilon parameter for the DBSCAN clustering algorithm.
#' @param minPts Integer. Minimum number of points to form a dense region in DBSCAN.
#' @param iter Number of iterations for the filtering process.
#' @param seeds Integer. Random seed for reproducibility.
#' @return A data frame containing the original group information with an added column
#'         `esg.prob`, which gives the probability of each group being identified as a
#'         extremely small group across all iterations.

#' @export

esg_filter <- function(group, n, r, adj_p=0.5, aprior=0.5, bprior=0.5, eps = 0.15, minPts = 4, iter = 1000, seeds = 07150321, ...) {
  set.seed(seeds)
  org.data <- data.frame(Group = group, n = n, r = r)

  outs <- lapply(1:iter, function(i) {
    data <- data.frame(subgroup = group, n = n, r = r)
    stop_loop <- FALSE
    small <- c()
    while (!stop_loop) {
      b.data <- data[is.na(match(data$subgroup, small)), ]
      adj_r <- purrr::pmap_vec(list(as.list(b.data$n), as.list(b.data$r)), function(n, r) adjust_responses(n, r, p = adj_p))
      jsd <- purrr::pmap_vec(list(as.list(b.data$n), as.list(b.data$r), as.list(adj_r)), function(x, y, z) jsd.two.beta(x, y, z, aprior = aprior, bprior = bprior))
      temp_data <- data.frame(subgroup = b.data$subgroup, n = b.data$n, r = b.data$r, adj_r = adj_r, jsd = jsd)
      temp_data$norm.n <- sqrt(temp_data$n / max(temp_data$n))
      temp_data$group <- dbscan::dbscan(temp_data[c("norm.n", "jsd")], eps = eps, minPts = minPts)$cluster
      small_group <- temp_data[temp_data$group == 0 & temp_data$n == min(temp_data$n) & temp_data$jsd == max(temp_data$jsd), "subgroup"]
      small <- c(small, small_group)
      a.data <- data[is.na(match(data$subgroup, small)), ]
      if (identical(b.data, a.data)) {
        stop_loop <- TRUE
      }
    }
    return(small)
  })

  df <- do.call(rbind, lapply(seq_along(outs), function(i) {
    data.frame(index = i, value = outs[[i]])
  }))

  results <- data.frame(table(df$value))
  names(results)[1] <- "Group"
  results$esg.prob <- results$Freq / iter
  final_results <- merge(org.data, results[, c(1, 3)], all = TRUE)
  final_results$esg.prob[is.na(final_results$esg.prob)] <- 0
  return(final_results)
}


#' Perform Gibbs Sampling for a Gaussian infinite mixture model
#'
#' This function performs Gibbs sampling on a Gaussian infinite mixture model
#' using a modified Chinese Restaurant Process (CRP) prior. It generates cluster
#' assignments for the data based on the Dirichlet Process.
#'
#' @param data A numeric vector (or n x 1 matrix) of data points to cluster.
#' @param lambda Numeric. Tuning parameters for the modified Chinese restaurant process mixed model.
#' @param alpha Numeric. The concentration parameter of the Dirichlet Process.
#' @param mu0 Numeric. Prior mean for the Gaussian distribution of cluster means.
#' @param sigma0 Numeric. Prior variance for the Gaussian distribution of cluster means.
#' @param sigma_y Numeric. Measurement error variance.
#' @param burnIn Integer. Number of burn-in iterations for Gibbs sampling.
#' @param maxIters Integer. The number of MCMC samples.
#'
#' @return A matrix where each row corresponds to the cluster assignments of the data points at a given iteration (after burn-in).
#'
#' @export

# Function to create Gibbs samples of Gaussian infinite mixture
GI.samp <- function(data,lambda=1,alpha=1e-20,mu0=0,sigma0=10,sigma_y=0.001,burnIn=1000,maxIters=1000) {

  tau0 <- solve(sigma0)
  tau_y <- solve(sigma_y)
  N <- length(data)
  data <- as.matrix(data)

  ##### Initialization
  # initialize the CRP Gibbs sampler
  z <- rep(1, N)             # initial cluster membership assignments
  n_k <- as.vector(table(z)) # initial data counts at each cluster
  Nclust <- length(n_k)      # initial number of clusters

  ##### CRP Gibbs sampler
  tables <- matrix(NA, nrow = N, ncol = burnIn+maxIters)
  for(iter in 1:(burnIn+maxIters)) {
    for(n in 1:N) {
      c_i <- z[n]
      n_k[c_i] <- n_k[c_i] - 1

      if( n_k[c_i] == 0 ) {
        n_k[c_i] <- n_k[Nclust]
        loc_z <- ( z == Nclust )
        z[loc_z] <- c_i
        n_k <- n_k[ -Nclust ]
        Nclust <- Nclust - 1
      }

      z[n] <- -1
      logp <- rep( NA, Nclust + 1 )
      for( c_i in 1:Nclust ) {
        tau_p <- tau0 + n_k[c_i] * tau_y
        sig_p <- solve(tau_p)
        loc_z <- which(z == c_i)
        sum_data <-  sum(data[z == c_i, ])
        mean_p <-  sig_p * (tau_y * sum_data + tau0 * mu0)
        logp[c_i] <- log(n_k[c_i]*lambda) + dnorm(data[n,], mean = mean_p, sd = sqrt(sig_p + sigma_y), log = TRUE) }

      logp[ Nclust+1 ] <- log(alpha) + dnorm(data[n,], mean = mu0, sd = sqrt(sigma0 + sigma_y), log = TRUE)
      max_logp <- max(logp)
      logp <- logp - max_logp
      loc_probs <- exp(logp)
      loc_probs <- loc_probs / sum(loc_probs)
      newz <- sample(1:(Nclust+1), 1, replace = TRUE, prob = loc_probs)
      if(newz == Nclust + 1) {
        n_k <- c(n_k, 0)
        Nclust <- Nclust + 1}
      z[n] <- newz
      n_k[newz] <- n_k[newz] + 1
    }
    tables[, iter] <- z
  }

  tables <- tables[,-(1:burnIn)]
  invisible(t(tables))

}


#' Perform Gibbs Sampling for a Finite Gaussian Mixture Model
#'
#' This function performs Gibbs sampling for a finite Gaussian mixture model.
#'
#' @param data A numeric vector or matrix representing the observed data points.
#' @param k Integer. Pre-specified the number of clusters.
#' @param mu0 Numeric. Prior mean for the Gaussian distribution of cluster means.
#' @param sigma0 Numeric. Prior variance for the Gaussian distribution of cluster means.
#' @param sigma_y Numeric. Measurement error variance.
#' @param a Numeric. Dirichlet prior parameters.
#' @param burnIn Integer. Number of burn-in iterations for Gibbs sampling.
#' @param maxIters Integer. The number of MCMC samples.
#'
#' @return A matrix where each row corresponds to the cluster assignments of the data points at a given iteration (after burn-in).
#'
#' @export

GF.samp <- function(data, k, burnIn=1000, maxIters=1000, mu0=0, sigma0=10, sigma_y=0.001, a=10000) {

  tau0 <- solve(sigma0)
  tau_y <- solve(sigma_y)
  N <- length(data)
  data <- as.matrix(data)

  ##### Initialization
  z <- rep(1:k, length.out = N)  # initial cluster membership assignments
  n_k <- as.vector(table(factor(z, levels = 1:k))) # initial data counts at each cluster

  ##### Gibbs sampler
  tables <- matrix(NA, nrow = N, ncol = burnIn+maxIters)
  for(iter in 1:(burnIn+maxIters)) {
    for(n in 1:N) {
      c_i <- z[n]
      n_k[c_i] <- n_k[c_i] - 1
      logp <- rep(NA, k)
      for( c_i in 1:k ) {
        tau_p <- tau0 + n_k[c_i] * tau_y
        sig_p <- solve(tau_p)
        loc_z <- which(z == c_i)
        sum_data <-  sum(data[loc_z, ])
        mean_p <-  sig_p * (tau_y * sum_data + tau0 * mu0)
        logp[c_i] <- log(n_k[c_i] + a/k) + dnorm(data[n,], mean = mean_p, sd = sqrt(sig_p + sigma_y), log = TRUE)
      }

      max_logp <- max(logp)
      logp <- logp - max_logp
      loc_probs <- exp(logp)
      loc_probs <- loc_probs / sum(loc_probs)
      newz <- sample(1:k, 1, replace = TRUE, prob = loc_probs)
      z[n] <- newz
      n_k[newz] <- n_k[newz] + 1
    }

    tables[, iter] <- z
  }

  tables <- tables[,-(1:burnIn)]
  invisible(t(tables))
}



# Function to calculate average silhouette index
avgsil <- function(data,samp){
  nsamp <- nrow(samp)
  avg <- vector("numeric", nsamp)
  diss <- as.matrix(dist(data))
  for (i in 1:nsamp) {
    sil <- cluster::silhouette(samp[i,], dmatrix=diss)
    if (!is.na(sil[1])) {avg[i] <- mean(sil[,3])}
    else {avg[i] <- NA}
  }
  if (!is.na(mean(avg, na.rm=T))) {return(mean(avg, na.rm=T))}
  else {return(0)}
}


# Function to calculate min and max of the number of clusters
minmax <- function(samp){
  out <- data.frame(min=min(apply(samp, 1, function(x) length(unique(x))),na.rm = T), max=max(apply(samp, 1, function(x) length(unique(x))),na.rm = T))
  return(out)
}

# Function to calculate Adjusted Rand Index
rand <- function(truecluster,samp){
  out <- data.frame(ARI=mean(apply(samp, 1, function(x) mclust::adjustedRandIndex(truecluster,x))))
  return(out)
}



#' Function to construct similarity matrix based on the response rates
#'
#' @param group A vector of subgroup labels.
#' @param n A vector representing the number of patients in each subgroup.
#' @param r A vector representing the number of responses in each subgroup.
#' @param kmin Integer. The minimum number of components (K) to consider in the Bayesian Finite Gaussian Mixture Model (BFGMM).
#' @param kmax Integer. The maximum number of components (K) to consider in the BFGMM.
#' @param lambda A numeric vector of tuning parameters for the modified Chinese restaurant process mixed model (m-CRPMM).
#' @param mu0 Numeric. The hyperprior mean for the BFGMM and m-CRPMM models.
#' @param sigma0 Numeric. The hyperprior variance for the BFGMM and m-CRPMM models.
#' @param sigma_y Numeric. Data measurement error for the BFGMM and m-CRPMM models.
#' @param alpha Numeric. The concentration parameter of the Dirichlet Process.
#' @param burnIn Integer. Number of burn-in iterations for Gibbs sampling.
#' @param maxIters Integer. The number of MCMC samples.
#'
#'
#' @return A list containing the following components:
#'   \item{data}{The original dataset of response rates provided to the function.}
#'   \item{summary}{A summary of the clustering analysis}
#'   \item{optmethod}{The optimal clustering method used.}
#'   \item{opttab}{MCMC samples from the optimal clustering method used.}
#'   \item{sm}{Similarity matrix based on clustering results.}
#'
#' @export

cluster.results <- function(group, n, r, alpha = 1, kmin=2, kmax=4, lambda=10^c(5,10,20,40,80), mu0=0, sigma0=10, sigma_y=0.001, burnIn=1000, maxIters=1000,...){
  res <- r/n
  y <- as.matrix(log((res+0.5/n)/(1-res+0.5/n)))  # empirical logits
  if (kmin>kmax) {kmin=kmax}
  kvec <- seq(kmin, kmax, by=1)
  lambda = lambda
  method <- c(paste0("BFGMM (K=", kvec,")"), "CRPMM (lambda=1)", paste0("m-CRPMM (lambda=", lambda,")"))
  gsamp.GF <- lapply(kvec, function(x) GF.samp(y, x, burnIn=burnIn, maxIters=maxIters))
  gsamp.GI <- lapply(c(1, lambda), function(x) GI.samp(y, x, burnIn=burnIn, maxIters=maxIters, alpha=alpha, mu0=mu0, sigma0=sigma0, sigma_y=sigma_y))
  gsamp <- c(gsamp.GF, gsamp.GI)

  names(gsamp) <- method
  sil <- data.frame(sil=purrr::map_vec(gsamp, function(x) avgsil(y,x)))
  max.sil <- sil |> dplyr::filter(sil==max(sil)) |> dplyr::slice(dplyr::n()) |> row.names()
  minmax <- purrr::map_df(gsamp, function(x) minmax(x))
  sum.tab <- dplyr::bind_cols(sil, minmax)
  table <- gsamp[[max.sil]]
  ngroup <- ncol(table)
  nruns <- nrow(table)
  sm <- matrix(0, ngroup, ngroup)
  for (i in 1:nruns)
  {
    rr <- table[i,]
    for (j in 1:ngroup)
    {
      for (k in 1:ngroup)
      {
        if (rr[j] == rr[k])
        {
          sm[j, k] <- sm[j, k] + 1
        }
      }
    }
  }
  sm.fin <- sm/nruns
  rownames(sm.fin) <- group
  colnames(sm.fin) <- group
  out <- list(data=data.frame(Group=group, n=n, r=r),
              summary=sum.tab,
              optmethod=max.sil,
              opttab=table,
              sm=sm.fin)
  return(out)

}


# Function to build JAG models
bmodel <- function() {
  mod1 <- "model{
  for (i in 1:numGroups) {
  y[i] ~ dbin(p[i],n[i])
  logit(p[i]) <- eta[i]
  rr[i] <- tau1 * w[i]
  eta[i] ~ dnorm(mu,rr[i])
  }
  #Priors
  mu ~ dnorm(mu0, tau2)
  tau1 ~ dgamma(alpha, beta)
  }"
  return(mod1)
}



# Function to construct similarity matrix based on the sample size
sm.d <- function(group, n) {
  ngroup <- length(n)
  d <- matrix(0, ngroup, ngroup)
  for (j in 1:ngroup) {
    for (i in 1:ngroup) {
      d[i, j] <- 1- max(0, (n[i] - n[j])/(n[j] + n[i]))
    }
  }
  rownames(d) <- group
  colnames(d) <- group
  return(d)
}
















