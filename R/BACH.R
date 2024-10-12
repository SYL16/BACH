#' BACH Function for Analyzing Concurrent Response Rate Data with Heterogeneity in Basket Trials
#'
#' @param group A vector of subgroup labels.
#' @param n A vector representing the number of patients in each subgroup.
#' @param r A vector representing the number of responses in each subgroup.
#' @param esfilter Logical. Whether to filter out extremely small subgroups. Default is TRUE.
#' @param escutoff Numeric. The cutoff probability used to identify extremely small subgroups.
#' @param kmin Integer. The minimum number of components (K) to consider in the Bayesian Finite Gaussian Mixture Model (BFGMM).
#' @param kmax Integer. The maximum number of components (K) to consider in the BFGMM.
#' @param lambda A numeric vector of tuning parameters for the modified Chinese restaurant process mixed model (m-CRPMM).
#' @param mu0 Numeric. The hyperprior mean for the BFGMM and m-CRPMM models.
#' @param sigma0 Numeric. The hyperprior variance for the BFGMM and m-CRPMM models.
#' @param sigma_y Numeric. Data measurement error for the BFGMM and m-CRPMM models.
#' @param w Numeric. A weight parameter (between 0 and 1) that determines the relative importance of subgroup similarity matrices \code{Cij} and \code{Dij}.
#' @param alpha Numeric. The concentration parameter (alpha) for the Dirichlet Process in m-CRPMM.
#' @param alpha1 Numeric. Prior shape parameter for borrowing strength in the Bayesian Hierarchical Model (BHM).
#' @param beta1 Numeric. Prior scale parameter for borrowing strength in the BHM.
#' @param tau2 Numeric. The hyperprior precision parameter for subgroup means in the BHM.
#' @param c.burnin Integer. The number of burn-in samples for the BFGMM and m-CRPMM models.
#' @param c.maxiters Integer. The number of MCMC iterations for the BFGMM and m-CRPMM models.
#' @param b.burnin Integer. The number of burn-in samples for the BHM model. Default is 5000.
#' @param b.maxiters Integer. The number of MCMC iterations for the BHM model. Default is 10000.
#' @param seeds Integer. Random seed for reproducibility.
#'
#' @return A list containing the following components:
#'   \item{optmethod}{The optimal clustering method used.}
#'   \item{mc}{Similarity matrix based on clustering results.}
#'   \item{mw}{Weighted similarity matrix combining clustering results and subgroup size similarities.}
#'   \item{Summary}{A data frame summarizing the posterior mean, 95\% credible interval, and other information for each subgroup.}
#'
#' @export


BACH <- function(group,
                 n,
                 r,
                 esfilter=T,
                 escutoff=0.8,
                 kmin=2,
                 kmax=4,
                 lambda=10^c(5,10,20,40,80),
                 mu0 = 0.2,
                 sigma0 = 10,
                 sigma_y=0.001,
                 w=0.95,
                 alpha = 1,
                 alpha1 = 50,
                 beta1 = 10,
                 tau2 = 0.1,
                 c.burnin=1000,
                 c.maxiters=1000,
                 b.burnin = 5000,
                 b.maxiters = 10000,
                 seeds=2024)
{

  numArm <- length(n)
  if (numArm != length(group) | numArm != length(r)) { stop("Length of group, n, and r are not equal.") }

  df.org <- data.frame(Group=group, No.Pat=n, No.Res=r, Obs.Rate=r/n)

  # Set random seed for reproducibility
  set.seed(seeds)

  # Apply extremely small filtering if esfilter is TRUE
  if (esfilter) {
    ddss_results <- esg_filter(group = group, n = n, r = r)
    filter_results <- ddss_results[ddss_results$esg.prob< escutoff, ]
    group <- filter_results$Group
    n <- filter_results$n
    r <- filter_results$r
  }

  group.R <- group
  n.R <- n
  r.R <- r
  ngroup <- length(n.R)

  c.results <-  cluster.results(group=group.R, n=n.R, r=r.R, seeds=seeds, kmin=kmin, kmax=kmax, lambda=lambda,
                                alpha = alpha, mu0=mu0, sigma0=sigma0, sigma_y=sigma_y,
                                burnIn=c.burnin, maxIters=c.maxiters)
  mc <- c.results$sm
  optmethod <- c.results$optmethod
  md <- sm.d(group=group.R, n=n.R)

  mw <- w*mc + (1-w)* md

  for (i in 1:ngroup)
  {
    cInd <- 1:ngroup
    t <- length(cInd)
    smP <- mw[i,]
    smP[smP<0.001] <- 0.001

    mydata <-
      list(
        y = r.R[cInd],
        n = n.R[cInd],
        w = smP,
        numGroups = length(cInd),
        mu0 = logit(mean(r.R[cInd] / n.R[cInd])),
        tau2 = tau2,
        alpha = alpha1,
        beta = beta1
      )

    mText1 <- bmodel()
    modelSpec1 <-textConnection(mText1)

    parameters <- c("p")

    jSeed <- floor(runif(1, 1, 10000))

    jags1 <- rjags::jags.model(
      modelSpec1,
      data = mydata,
      n.chains = 4,
      n.adapt = b.maxiters / 4,
      quiet = TRUE,
      inits = list(.RNG.name = "base::Wichmann-Hill",
                   .RNG.seed = jSeed)
    )

    # Burn-in period
    update(jags1, n.iter = b.burnin, progress.bar	= "none")

    MCRes1 <- rjags::coda.samples(jags1,
                           parameters,
                           n.iter = b.maxiters,
                           verbose = FALSE,
                           progress.bar = "none",
                           thin = 1)

    samples <- MCRes1[[1]]
    sampledP <- samples[, i]

    if (i == 1)
    {
      allPost <- sampledP
    } else{
      allPost <- cbind(allPost, sampledP)
    }

  }

  colnames(allPost) <- group.R

  post.summ <- data.frame(allPost) |>
    tidyr::pivot_longer(everything(),names_to = "Group", values_to = "p") |>
    dplyr::group_by(Group) |>
    dplyr::summarise(
      Post.Mean=round(mean(p),3),
      Lower.Bound=round(quantile(p, 0.025),3),
      Upper.Bound=round(quantile(p, 0.975),3))

  outs <- suppressMessages(df.org |>
                             dplyr::left_join(post.summ) |>
                             dplyr::mutate(Comments=ifelse(Group %in% group, NA, "Identified as an extremely small group.")) |>
                             dplyr::arrange(No.Pat))

  if (esfilter) { results <- list(summary_ds=ddss_results,
                                  summary_cluster=c.results$summary,
                                  optmethod=optmethod, opttab=c.results$opttab,
                                  mc = mc, mw = mw,
                                  post_summary = outs, post=allPost)}
  else { results <- list(summary_cluster=c.results$summary,
                         optmethod=optmethod, opttab=c.results$opttab,
                         mc = mc, mw = mw,
                         summary_post = outs, post=allPost) }

  return(results)

}



