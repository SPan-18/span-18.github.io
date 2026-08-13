rm(list = ls())

## Section 2.1: Conjugate Bayesian Gaussian spatial model
library("spStack")
data(simGaussian)

# setup prior list
muBeta <- c(0, 0)
VBeta <- diag(1E3, 2)
sigmaSqIGa <- 2
sigmaSqIGb <- 2
prior_list <- list(beta.norm = list(muBeta, VBeta),
                   sigma.sq.ig = c(sigmaSqIGa, sigmaSqIGb))

# exact sampling from the posterior
set.seed(1729)
mod1 <- spLMexact(y ~ x1, data = simGaussian,
                  coords = as.matrix(simGaussian[, c("s1", "s2")]),
                  cor.fn = "matern",
                  priors = prior_list,
                  spParams = list(phi = 3, nu = 0.75),
                  noise_sp_ratio = 0.8, 
                  n.samples = 1000, verbose = TRUE)

beta.post <- mod1$samples$beta
rownames(beta.post) <- mod1[["X.names"]]
print(t(apply(beta.post, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))))

## Section 2.2: Conjugate Bayesian non-Gaussian spatial model
data(simPoisson)
mod2 <- spGLMexact(y ~ x1, data = simPoisson, family = "poisson",
                   coords = as.matrix(simPoisson[, c("s1", "s2")]),
                   cor.fn = "matern", spParams = list(phi = 4, nu = 0.4),
                   boundary = 0.5, n.samples = 1000, verbose = TRUE)
beta.post <- mod2$samples$beta
rownames(beta.post) <- mod1[["X.names"]]
print(t(apply(beta.post, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))))

## Section 5: Illustrations

# test exact vs PSIS loopd
set.seed(1729)
mod3 <- spLMexact(y ~ x1, data = simGaussian,
                  coords = as.matrix(simGaussian[, c("s1", "s2")]),
                  cor.fn = "matern", spParams = list(phi = 3, nu = 0.75),
                  noise_sp_ratio = 0.8, n.samples = 100,
                  loopd = TRUE, loopd.method = "exact", verbose = FALSE)

mod4 <- spLMexact(y ~ x1, data = simGaussian,
                  coords = as.matrix(simGaussian[, c("s1", "s2")]),
                  cor.fn = "matern", spParams = list(phi = 3, nu = 0.75),
                  noise_sp_ratio = 0.8, n.samples = 100,
                  loopd = TRUE, loopd.method = "PSIS", verbose = FALSE)

loopd_exact <- mod3$loopd
loopd_psis <- mod4$loopd
loopd_df <- data.frame(exact = loopd_exact, psis = loopd_psis)
library(ggplot2)
ggplot(data = loopd_df, aes(x = exact)) +
  geom_point(aes(y = psis), size = 1, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", alpha = 0.5) +
  xlab("Exact") + ylab("PSIS") + theme_bw() +
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(),
        aspect.ratio = 1)

# library(microbenchmark)
# et1 <- microbenchmark(
#   spLMexact(y ~ x1, data = simGaussian,
#             coords = as.matrix(simGaussian[, c("s1", "s2")]),
#             cor.fn = "matern", spParams = list(phi = 3, nu = 0.75),
#             noise_sp_ratio = 0.8, n.samples = 100,
#             loopd = TRUE, loopd.method = "exact", verbose = FALSE),
#   times = 10
# )
# 
# et2 <- microbenchmark(
#   spLMexact(y ~ x1, data = simGaussian,
#             coords = as.matrix(simGaussian[, c("s1", "s2")]),
#             cor.fn = "matern", spParams = list(phi = 3, nu = 0.75),
#             noise_sp_ratio = 0.8, n.samples = 100,
#             loopd = TRUE, loopd.method = "PSIS", verbose = FALSE),
#   times = 10
# )
# 
# median(et1$time) / median(et2$time)

# stacking for Gaussian data
mod5 <- spLMstack(y ~ x1, data = simGaussian,
                  coords = as.matrix(simGaussian[, c("s1", "s2")]),
                  cor.fn = "matern",
                  params.list = list(phi = c(1.5, 3, 5),
                                     nu = c(0.5, 1.5),
                                     noise_sp_ratio = c(0.5, 1.5)),
                  n.samples = 1000, loopd.method = "PSIS",
                  parallel = FALSE, solver = "ECOS", verbose = TRUE)

print(mod5[["solver.status"]])

# sample from the stacked posterior
post_samps <- stackedSampler(mod5)
post_beta <- post_samps[["beta"]]
summary_beta <- t(apply(post_beta, 1, 
                        function(x) quantile(x, c(0.025, 0.5, 0.975))))
rownames(summary_beta) <- mod5[["X.names"]]
print(summary_beta)

post_z <- post_samps[["z"]]
post_z_summ <- t(apply(post_z, 1, 
                       function(x) quantile(x, c(0.025, 0.5, 0.975))))
z_combn <- data.frame(z = simGaussian[["z_true"]], zL = post_z_summ[, 1],
                      zM = post_z_summ[, 2], zU = post_z_summ[, 3])

errbar_z <- ggplot(data = z_combn, aes(x = z)) +
  geom_point(aes(y = zM), size = 0.5, color = "darkblue", alpha = 0.5) +
  geom_errorbar(aes(ymin = zL, ymax = zU), width = 0.05, alpha = 0.15, 
                color = "skyblue") + 
  geom_abline(slope = 1, intercept = 0, color = "red") +
  xlab("True z") + ylab("Stacked posterior of z") + theme_bw() +
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        aspect.ratio = 1)
errbar_z

postmedian_z <- apply(post_z, 1, median)
simGaussian$z_hat <- postmedian_z
plot_z <- surfaceplot2(simGaussian, coords_name = c("s1", "s2"),
                       var1_name = "z_true", var2_name = "z_hat")
plot_z[[1]] <- plot_z[[1]] + xlab(bquote(s[1])) + ylab(bquote(s[2]))
plot_z[[2]] <- plot_z[[2]] + xlab(bquote(s[1])) + ylab(bquote(s[2]))
library(ggpubr)
ggarrange(plotlist = plot_z, common.legend = TRUE, legend = "right")

## Section 5.2: Spatial Poisson count data
data(simPoisson)

pois_raw <- ggplot(simPoisson, aes(x = s1, y = s2)) +
  geom_point(aes(color = y, alpha = log(y)), size = 1.5) +
  scale_color_distiller(palette = "RdYlGn", direction = -1, 
                        label = function(x) sprintf("%.0f", x)) +
  guides(alpha = 'none') + theme_bw() + xlab(bquote(s[1])) + ylab(bquote(s[2])) +
  theme(axis.ticks = element_line(linewidth = 0.25),
        panel.background = element_blank(), panel.grid = element_blank(),
        legend.title = element_text(size = 10, hjust = 0.25),
        legend.box.just = "center", aspect.ratio = 1)

cand_list <- list(phi = c(3, 7, 10), nu = c(0.5, 1.5), boundary = c(0.5, 0.6))
loopd_settings <- list(method = "CV", CV.K = 10, nMC = 1000)

library("future")
plan("multicore", workers = 6)
set.seed(1729)
mod6 <- spGLMstack(y ~ x1, data = simPoisson, family = "poisson",
                   coords = as.matrix(simPoisson[, c("s1", "s2")]), 
                   cor.fn = "matern", params.list = cand_list, 
                   n.samples = 1000, loopd.controls = loopd_settings,
                   parallel = TRUE, solver = "ECOS", verbose = TRUE)
plan("sequential")

post_samps <- stackedSampler(mod6)
post_beta <- post_samps[["beta"]]
summary_beta <- t(apply(post_beta, 1, 
                        function(x) quantile(x, c(0.025, 0.5, 0.975))))
rownames(summary_beta) <- mod6[["X.names"]]
print(summary_beta)

post_z <- post_samps[["z"]]
post_z_summ <- t(apply(post_z, 1, 
                       function(x) quantile(x, c(0.025, 0.5, 0.975))))
z_combn <- data.frame(z = simPoisson[["z_true"]], zL = post_z_summ[, 1],
                      zM = post_z_summ[, 2], zU = post_z_summ[, 3])

ggplot(data = z_combn, aes(x = z)) +
  geom_point(aes(y = zM), size = 0.5, color = "darkblue", alpha = 0.5) +
  geom_errorbar(aes(ymin = zL, ymax = zU), width = 0.05, alpha = 0.15, 
                color = "skyblue") + 
  geom_abline(slope = 1, intercept = 0, color = "red") +
  xlab("True z") + ylab("Stacked posterior of z") + theme_bw() +
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        aspect.ratio = 1)

postmedian_z <- apply(post_z, 1, median)
simPoisson$z_hat <- postmedian_z
plot_z2 <- surfaceplot2(simPoisson, coords_name = c("s1", "s2"),
                        var1_name = "z_true", var2_name = "z_hat")
plot_z2[[1]] <- plot_z2[[1]] + xlab(bquote(s[1])) + ylab(bquote(s[2]))
plot_z2[[2]] <- plot_z2[[2]] + xlab(bquote(s[1])) + ylab(bquote(s[2]))
plot_z2_comb <- ggarrange(plotlist = plot_z2, common.legend = TRUE, 
                          legend = "right")

# Figure 5
ggarrange(pois_raw, plot_z2_comb, labels = c("(a)", "(b)"),
          widths = c(1, 1.8))

## Section 5.2: Spatial Binomial count data
data(simBinom)

library("future")
plan("multicore", workers = 6)
cand_list <- list(phi = c(2, 5, 8), nu = c(1.0, 1.5), boundary = c(0.5, 0.6))
mod7 <- spGLMstack(cbind(y, n_trials) ~ x1, data = simBinom, 
                   family = "binomial",
                   coords = as.matrix(simBinom[, c("s1", "s2")]), 
                   cor.fn = "matern", params.list = cand_list, 
                   n.samples = 1000, loopd.controls = loopd_settings,
                   parallel = TRUE, solver = "ECOS", verbose = TRUE)
plan("sequential")

post_samps <- stackedSampler(mod7)
post_beta <- post_samps[["beta"]]
summary_beta <- t(apply(post_beta, 1, 
                        function(x) quantile(x, c(0.025, 0.5, 0.975))))
rownames(summary_beta) <- mod7[["X.names"]]
print(summary_beta)

post_z <- post_samps[["z"]]
post_z_summ <- t(apply(post_z, 1, 
                       function(x) quantile(x, c(0.025, 0.5, 0.975))))
z_combn <- data.frame(z = simBinom[["z_true"]], zL = post_z_summ[, 1],
                      zM = post_z_summ[, 2], zU = post_z_summ[, 3])

ggplot(data = z_combn, aes(x = z)) +
  geom_point(aes(y = zM), size = 0.5, color = "darkblue", alpha = 0.5) +
  geom_errorbar(aes(ymin = zL, ymax = zU), width = 0.05, alpha = 0.15, 
                color = "skyblue") + 
  geom_abline(slope = 1, intercept = 0, color = "red") +
  xlab("True z") + ylab("Stacked posterior of z") + theme_bw() +
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        aspect.ratio = 1)

postmedian_z <- apply(post_z, 1, median)
simBinom$z_hat <- postmedian_z
plot_z3 <- surfaceplot2(simBinom, coords_name = c("s1", "s2"),
                        var1_name = "z_true", var2_name = "z_hat")
plot_z3[[1]] <- plot_z3[[1]] + xlab(bquote(s[1])) + ylab(bquote(s[2]))
plot_z3[[2]] <- plot_z3[[2]] + xlab(bquote(s[1])) + ylab(bquote(s[2]))
plot_z3_comb <- ggarrange(plotlist = plot_z3, common.legend = TRUE, 
                          legend = "right")
plot_z3_comb

## Section 5.2: Spatial binary data
data(simBinary)

binary_raw <- ggplot(simBinary, aes(x = s1, y = s2)) +
  geom_point(aes(color = factor(y)), alpha = 0.75, , size = 1.5) +
  scale_color_manual(values = c("red", "blue"), labels = c("0", "1")) +
  guides(alpha = 'none') + theme_bw() + labs(color = "y") +
  theme(axis.ticks = element_line(linewidth = 0.25),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = 10, hjust = 0.25),
        legend.box.just = "center", aspect.ratio = 1)

library("future")
plan("multicore", workers = 6)
cand_list <- list(phi = c(2, 5, 8), nu = c(0.75, 1.25), boundary = c(0.5, 0.6))
mod8 <- spGLMstack(y ~ x1, data = simBinary, family = "binary",
                   coords = as.matrix(simBinary[, c("s1", "s2")]), 
                   cor.fn = "matern", params.list = cand_list, 
                   n.samples = 1000, loopd.controls = loopd_settings,
                   parallel = TRUE, solver = "ECOS", verbose = TRUE)
plan("sequential")

post_samps <- stackedSampler(mod8)
post_beta <- post_samps[["beta"]]
summary_beta <- t(apply(post_beta, 1, 
                        function(x) quantile(x, c(0.025, 0.5, 0.975))))
rownames(summary_beta) <- mod8[["X.names"]]
print(summary_beta)

post_z <- post_samps[["z"]]
post_z_summ <- t(apply(post_z, 1, 
                       function(x) quantile(x, c(0.025, 0.5, 0.975))))
z_combn <- data.frame(z = simBinary[["z_true"]], zL = post_z_summ[, 1],
                      zM = post_z_summ[, 2], zU = post_z_summ[, 3])

ggplot(data = z_combn, aes(x = z)) +
  geom_point(aes(y = zM), size = 0.5, color = "darkblue", alpha = 0.5) +
  geom_errorbar(aes(ymin = zL, ymax = zU), width = 0.05, alpha = 0.15, 
                color = "skyblue") + 
  geom_abline(slope = 1, intercept = 0, color = "red") +
  xlab("True z") + ylab("Stacked posterior of z") + theme_bw() +
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        aspect.ratio = 1)

## Appendix B
set.seed(1729)
tol <- 1E-12
n <- 10
A <- matrix(rnorm(n^2), n, n)
A <- crossprod(A)
cholA <- chol(A)

# test rank-1 update
v <- 1:n
APlusvvT <- A + tcrossprod(v)
cholA1 <- t(chol(APlusvvT))
cholA2 <- cholUpdateRankOne(cholA, v, alpha = 1, beta = 1, lower = F)
print(max(abs(cholA1 - cholA2)) < tol)

# test single row-deletion update
ind <- 2
A1 <- A[-ind, -ind]
cholA3 <- t(chol(A1))
cholA4 <- cholUpdateDel(cholA, del.index = ind, lower = F)
print(max(abs(cholA3 - cholA4)) < tol)

# test block-deletion update
start_ind <- 2
end_ind <- 6
del_ind <- c(start_ind:end_ind)
A2 <- A[-del_ind, -del_ind]
cholA5 <- t(chol(A2))
cholA6 <- cholUpdateDelBlock(cholA, del.start = start_ind, 
                             del.end = end_ind, lower = F)
print(max(abs(cholA5 - cholA6)) < tol)

sessionInfo()
