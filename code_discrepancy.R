## ----setup, include=FALSE--------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, dev = "tikz", cache = TRUE)


## ---- results="hide", message=FALSE, warning=FALSE, cache=FALSE------------------

# PRELIMINARY ------------------------------------------------------------------

# Function to read in all required packages in one go:
loadPackages <- function(x) {
  for(i in x) {
    if(!require(i, character.only = TRUE)) {
      install.packages(i, dependencies = TRUE)
      library(i, character.only = TRUE)
    }
  }
}

theme_AP <- function() {
  theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "transparent",
                                           color = NA),
          legend.margin=margin(0, 0, 0, 0),
          legend.box.margin=margin(-7,-7,-7,-7), 
          legend.key = element_rect(fill = "transparent",
                                    color = NA), 
          strip.background = element_rect(fill = "white"))
}

# Load the packages
loadPackages(c("sensobol", "data.table", "tidyverse", "parallel", 
               "RcppAlgos", "scales", "doParallel", "benchmarkme", 
               "cowplot"))


## ----savage_fun------------------------------------------------------------------

# SAVAGE SCORES FUNCTION -------------------------------------------------------

savage_scores <- function(x) {
  true.ranks <- rank(-x)
  p <- sort(1 / true.ranks)
  mat <- matrix(rep(p, length(p)), nrow = length(p), byrow = TRUE)
  mat[upper.tri(mat)] <- 0
  out <- sort(rowSums(mat), decreasing = TRUE)[true.ranks]
  return(out)
}


## ----discrepancy_fun, dependson="savage_scores"----------------------------------

# DISCREPANCY FUNCTION ---------------------------------------------------------

# This is based on the function discrepancyCriteria_cpp of
# the sensitivity package

# Function to rescale ---------------------------
rescale_fun <- function(x) (x - min(x)) / (max(x) - min(x))

# Source cpp code -------------------------------

cpp_functions <- c("cpp_functions.cpp", "L2star_functions.cpp", 
                   "L2_functions.cpp", "L2centered_functions.cpp", 
                   "L2wraparound_functions.cpp")

for(i in 1:length(cpp_functions)) {
  Rcpp::sourceCpp(cpp_functions[i])
}

# Function --------------------------------------
discrepancy_fun <- function (design, type) {

  X <- as.matrix(design)
  dimension <- ncol(X)
  n <- nrow(X)
  
  # Reescale if needed-------
  
  if (min(X) < 0 || max(X) > 1) {
    
    X <- apply(X, 2, rescale_fun)
  }
  
  # Compute discrepancy
  
  if (type == "symmetric") {
    
    P <- 1 + 2 * X - 2 * X^2
    s1 <- DisS2_Rowprod(t(P), dimension)
    s2 <- DisS2_Crossprod(c(t(X)), dimension)
    R <- sqrt(((4/3)^dimension) - ((2/n) * s1) + ((2^dimension/n^2) * s2))
  
  } else if (type == "star") {
    
    dL2 <- DisL2star_Crossprod(t(X), dimension)
    R <- sqrt(3^(-dimension) + dL2)
    
  } else if (type == "L2") {
    
    P <- X * (1 - X)
    s1 <- DisL2_Rowprod(t(P), dimension)
    s2 <- DisL2_Crossprod(c(t(X)), dimension)
    R <- sqrt(12^(-dimension) - (((2^(1 - dimension))/n) * s1) + ((1/n^2) * s2))
    
  } else if (type == "centered") {
    
    P <- 1 + 0.5 * abs(X - 0.5) - 0.5 * (abs(X - 0.5)^2)
    s1 <- DisC2_Rowprod(t(P), dimension)
    s2 <- DisC2_Crossprod(c(t(X)), dimension)
    R <- sqrt(((13/12)^dimension) - ((2/n) * s1) + ((1/n^2) * s2))
    
  } else if (type == "wraparound") {
    
    s1 <- DisW2_Crossprod(t(X), dimension)
    R <- sqrt(-(4/3)^dimension + (1/n^2) * s1)
    
  }
  return(R)
}


# Discrepancy function --------------------------
discrepancy <- function(mat, y, params, type) {
  value <- sapply(1:ncol(mat), function(j) {
    design <- cbind(mat[, j], y) 
    value <- discrepancy_fun(design = design, type = type)
  })
  return(value)
}


## ----jansen_fun------------------------------------------------------------------

# FUNCTION TO COMPUTE JANSEN TI ------------------------------------------------

jansen_ti <- function(d, N, params) {
  m <- matrix(d, nrow = N)
  k <- length(params)
  Y_A <- m[, 1]
  Y_AB <- m[, -1]
  f0 <- (1 / length(Y_A)) * sum(Y_A)
  VY <- 1 / length(Y_A) * sum((Y_A - f0) ^ 2)
  value <- (1 / (2 * N) * Rfast::colsums((Y_A - Y_AB) ^ 2)) / VY
  return(value)
}


## ----replicas--------------------------------------------------------------------

# FUNCTION TO CREATE REPLICAS OF SAMPLE MATRIX ---------------------------------

# For discrepancy
CutBySize <- function(m, block.size, nb = ceiling(m / block.size)) {
  int <- m / nb
  upper <- round(1:nb * int)
  lower <- c(1, upper[-nb] + 1)
  size <- c(upper[1], diff(upper))
  cbind(lower, upper, size)
}

sobol_replicas <- function(matrices, N, params, replicas,...) {
  mat <- sobol_matrices(matrices = matrices, N = N * replicas, 
                        params = params,...)
  index <- CutBySize(nrow(mat), block.size = N)
  out <- list()
  for (i in 1:nrow(index)) {
    out[[i]] <- mat[index[i, "lower"]:index[i, "upper"], ]
  }
  return(out)
}

# For Sobol' indices
scrambled_sobol <- function(A, B) {
  X <- rbind(A, B)
  for(i in 1:ncol(A)) {
    AB <- A
    AB[, i] <- B[, i]
    X <- rbind(X, AB)
  }
  AB <- rbind(A, X[((2*nrow(A)) + 1):nrow(X), ])
  return(AB)
}
# Function to create replicas of the A, B and AB matrices
sobol_replicas_A_AB <- function(N, params, replicas) {
  k <- length(params)
  df <- randtoolbox::sobol(n = N * replicas, dim = k * 2)
  indices <- CutBySize(nrow(df), nb = replicas)
  X <- A <- B <- out <- list()
  for(i in 1:nrow(indices)) {
    lower <- indices[i, "lower"]
    upper <- indices[i, "upper"]
    X[[i]] <- df[lower:upper, ]
  }
  for(i in seq_along(X)) {
    A[[i]] <- X[[i]][, 1:k]
    B[[i]] <- X[[i]][, (k + 1) : (k * 2)]
  }
  for(i in seq_along(A)) {
    out[[i]] <- scrambled_sobol(A[[i]], B[[i]])
  }
  return(out)
}

## ----correlation, dependson=c("replicas", "jansen_fun", "discrepancy_fun", "savage_fun")----

# FUNCTION TO CHECK CORRELAITON BETWEEN DISCREPANCIES AND SAVAGE SCORES --------

cor_fun <- function(N = N, params = params, replicas = replicas, 
                    type = "symmetric", approach = NULL,
                    model_fun = NULL, true_scores = NULL,...) {
  
  if (is.null(approach)) {
    stop("approach should be either discrepancy or sobol")
  }
  
  if (is.null(model_fun) | is.null(true_scores)) {
    stop("A function and the true ranks should be provided")
  }
  
  if (approach == "discrepancy") {
    
    mat <- sobol_replicas(matrices = "A", N = N, params = params, 
                          replicas = replicas,...)
    y <- lapply(mat, model_fun)
    
    out <- list()
    for(i in 1:length(y)) {
      out[[i]] <- discrepancy(mat = mat[[i]], y = y[[i]], params = params, 
                              type = type) %>%
        savage_scores()
    }
    
  } else if (approach == "sobol") {
    
    mat <- sobol_replicas_A_AB(N = N, params = params, replicas = replicas)
    y <- lapply(mat, model_fun)
    out <- lapply(y, function(x) jansen_ti(d = x, N = N, params = params) %>%
                    savage_scores()) 
  }
  
  final <- unlist(lapply(out, function(x) cor(x, true_scores)))
  
  return(final)
}


## ----bratley_check---------------------------------------------------------------

# CHECK SOBOL' INDICES AND SCATTERPLOT FOR SEVERAL FUNCTIONS ------------------

N <- 2^14 # Base sample size
matrices <- c("A", "B", "AB")
R <- 10^3 # Bootstrap replicas

models <- c("Bratley et al. 1992", "Oakland and O'Hagan 2004")
fun_list <- list(bratley1992_Fun, oakley_Fun)
names(fun_list) <- models

out <- ind <- savage.rank <- list()
for (i in models) {
  
  if (i == "Bratley et al. 1992") {
    params <- paste("X", 1:6, sep = "")
    
  } else if (i == "Oakland and O'Hagan 2004") {
    params <- paste("X", 1:15, sep = "")
  }
  
  mat <- sobol_matrices(matrices = matrices, N = N, params = params)
  y <- fun_list[[i]](mat)
  ind[[i]] <- sobol_indices(matrices = matrices, Y = y, N = N, params = params, 
                       boot = TRUE, R = R)
  savage.rank[[i]] <- ind[[i]]$results[sensitivity == "Ti"] %>%
    .[, savage_scores(original)]
  out[[i]] <-plot_scatter(data = mat, N = N, Y = y, params = params) 
}

lapply(ind, plot)


## ----model, dependson="correlation"----------------------------------------------

# THE MODEL -------------------------------------------------------------------

matrices <- "A"
replicas <- 150 # Number of replicas of the sample matrix

output <-  list()

for(i in models) {
  
  if (i == "Bratley et al. 1992") {
    params <- paste("X", 1:6, sep = "")
    
  } else if (i == "Oakland and O'Hagan 2004") {
    params <- paste("X", 1:15, sep = "")

  }
  
  for (j in c("discrepancy", "sobol")) {
    
    if (j == "discrepancy") {
      
      if(i == "Bratley et al. 1992") {
        sample.sizes <- seq(4, 150, by = 4)
        final.sample.sizes <- sample.sizes
      } else {
        sample.sizes <- seq(4, 320, by = 4)
        final.sample.sizes <- sample.sizes
      }

    } else if (j == "sobol") {
      sample.sizes <- 2:20
      final.sample.sizes <- sample.sizes * (length(params) + 1)
      
    }
    output[[i]][[j]] <- mclapply(sample.sizes, function(x) 
      cor_fun(N = x, params = params, replicas = replicas, model_fun = fun_list[[i]], 
              true_scores = savage.rank[[i]], approach = j, type = "symmetric",
              scrambling = 1), 
      mc.cores = detectCores() * 0.75)
    
    names(output[[i]][[j]]) <- final.sample.sizes

  }
}

# Arrange output -------------------
results <- lapply(output, function(x) lapply(x, function(y) 
  lapply(y, function(z) data.table(z))))  %>%
  lapply(., function(x) lapply(x, function(y) rbindlist(y, idcol = "Model.runs"))) %>%
  lapply(., function(x) rbindlist(x, idcol = "Approach")) %>%
  rbindlist(., idcol = "Function") %>%
  .[, Model.runs:= as.numeric(Model.runs)] %>%
  .[, Approach:= ifelse(Approach == "discrepancy", "Discrepancy", "Jansen estimator")] %>%
  na.omit()



## ----plot, dependson="model", warning=FALSE, fig.height=4, fig.cap="Correlation between the Savage scores obtained with $T_i$ (Jansen estimator) and those obtained with discrepancy."----

# Plot the results ------------------------------------------------------------

# New facet label names for supp variable
supp.labs <- c("Symmetric $L2$ discrepancy", "Jansen estimator")
names(supp.labs) <- c("Discrepancy", "Jansen estimator")

ggplot(results, aes(Model.runs, z, group = Model.runs)) +
  geom_boxplot(position = "identity", outlier.size = 0.5) +
  labs(x = "Nº of model runs", y = "Correlation") +
  facet_grid(Approach~Function, 
             scale = "free_x", space = "free", 
             labeller = labeller(Approach = supp.labs)) +
  theme_AP()

## ----metamodel_fun---------------------------------------------------------------

# RANDOM FUNCTIONS AND DISTRIBUTIONS ------------------------------------------

# Functions ------------------------------------
function_list <- list(
  Linear = function(x) x,
  Quadratic = function(x) x ^ 2,
  Cubic = function(x) x ^ 3,
  Exponential = function(x) exp(1) ^ x / (exp(1) - 1),
  Periodic = function(x) sin(2 * pi * x) / 2,
  Discontinuous = function(x) ifelse(x > 0.5, 1, 0),
  Non.monotonic = function(x) 4 * (x - 0.5) ^ 2,
  Inverse = function(x) (10 - 1 / 1.1) ^ -1 * (x + 0.1) ^ - 1, 
  No.effect = function(x) x * 0, 
  Trigonometric = function(x) cos(x), 
  Piecewise.large = function(x) ((-1) ^ as.integer(4 * x) * 
                                   (0.125 - (x %% 0.25)) + 0.125),
  Piecewise.small = function(x) ((-1) ^ as.integer(32 * x) * 
                                   (0.03125 - 2 * (x %% 0.03125)) + 0.03125) / 2,
  Oscillation = function(x) x ^ 2 - 0.2 * cos(7 * pi * x)
)

# Random distributions  -------------------------
sample_distributions <- list(
  "uniform" = function(x) x,
  "normal" = function(x) qnorm(x, 0.5, 0.15),
  "beta" = function(x) qbeta(x, 8, 2),
  "beta2" = function(x) qbeta(x, 2, 8),
  "beta3" = function(x) qbeta(x, 2, 0.8),
  "beta4" = function(x) qbeta(x, 0.8, 2),
  "logitnormal" = function(x) logitnorm::qlogitnorm(x, 0, 3.16)
  # Logit-normal, Bates too?
)

random_distributions <- function(X, phi) {
  names_ff <- names(sample_distributions)
  if(!phi == length(names_ff) + 1) {
    out <- sample_distributions[[names_ff[phi]]](X)
  } else {
    temp <- sample(names_ff, ncol(X), replace = TRUE)
    out <- sapply(seq_along(temp), function(x) 
      sample_distributions[[temp[x]]](X[, x]))
  }
  return(out)
}

## ----model_fun, dependson=c("metamodel_fun", "savage_fun", "discrepancy_fun", "jansen_fun")----

# DEFINE MODEL ----------------------------------------------------------------

model_fun <- function(tau, epsilon, base.sample.size, cost.discrepancy, phi, k) {
  
  params <- paste("X", 1:k, sep = "")
  
  if (tau == 1) {
    
    type <- "R"
    
  } else if (tau == 2) {
    
    type <- "QRN"
  }
  set.seed(epsilon)
  discrepancy.mat <- sobol_matrices(matrices = "A", N = cost.discrepancy, 
                                    params = params, 
                                    type = type)
  
  set.seed(epsilon)
  jansen.mat <- sobol_matrices(matrices = c("A", "AB"), N = base.sample.size, 
                               params = params, 
                               type = type)
  
  all.matrices <- rbind(discrepancy.mat, jansen.mat)
  set.seed(epsilon)
  transformed.all.matrices <- random_distributions(X = all.matrices, phi = phi)
  
  y <- sensobol::metafunction(data = transformed.all.matrices, epsilon = epsilon)
  
  type <- c("symmetric", "star", "L2", "centered", "wraparound")
  
  discrepancy.value <- lapply(type, function(x) 
    discrepancy(mat = all.matrices[1:cost.discrepancy, ], 
                y = y[1:cost.discrepancy], params = params, type = x))
    
  
  jansen.value <- jansen_ti(d = y[(cost.discrepancy + 1):length(y)],
                            N = base.sample.size, params = params)
  
  savage.discrepancy <- lapply(discrepancy.value, savage_scores)
  jansen.discrepancy <- savage_scores(jansen.value)
  
  out <- lapply(savage.discrepancy, function(x) 
    cor(x, jansen.discrepancy))
    
  return(out)
}


## ----matrix----------------------------------------------------------------------

# CREATE SAMPLE MATRIX -------------------------------------------------------

N <- 2^2
params <- c("epsilon", "phi", "k", "tau", "base.sample.size")
mat <- sobol_matrices(matrices = "A", N = N, params = params)

# Define distributions --------------------------

mat[, "epsilon"] <- floor(qunif(mat[, "epsilon"], 1, 200))
mat[, "phi"] <- floor(mat[, "phi"] * 8) + 1
mat[, "k"] <- floor(qunif(mat[, "k"], 3, 50))
mat[, "tau"] <- floor(mat[, "tau"] * 2) + 1
mat[, "base.sample.size"] <- floor(qunif(mat[, "base.sample.size"], 10, 100))

cost.jansen <- mat[, "base.sample.size"] * (mat[, "k"] + 1)
cost.discrepancy <- cost.jansen

final.mat <- cbind(mat, cost.jansen, cost.discrepancy)


## ----run_model_fun, dependson="model_fun"----------------------------------------

# RUN MODEL  ------------------------------------------------------------------

y <- mclapply(1:nrow(final.mat), function(i) {
  model_fun(tau = final.mat[i, "tau"],
            epsilon = final.mat[i, "epsilon"], 
            base.sample.size = final.mat[i, "base.sample.size"], 
            cost.discrepancy = final.mat[i, "cost.discrepancy"], 
            phi = final.mat[i, "phi"], 
            k = final.mat[i, "k"])}, 
  mc.cores = floor(detectCores() * 0.75))

# ARRANGE DATA -----------------------------------------------------------------

y <- lapply(y, unlist)
output <- data.table(do.call(rbind, y))
discrepancy_methods <- c("symmetric", "star", "L2", "centered", "wraparound")
colnames(output) <- discrepancy_methods

final.output <- data.table(cbind(final.mat, output))

## ----plot_uncertainty, dependson="run_model_fun", warning = FALSE, fig.height=3, fig.width=3----

# PLOT UNCERTAINTY ------------------------------------------------------------

# New facet label names for supp variable
supp.labs <- c("Symmetric", "Star", "$L_2$", "Centered", "Wrap-around")
names(supp.labs) <- discrepancy_methods

a <- melt(final.output, measure.vars = discrepancy_methods) %>%
  .[, variable:= factor(variable, levels = c("star", "L2", "centered", 
                                                "wraparound", "symmetric"))] %>%
  ggplot(., aes(cost.discrepancy, k, color = value)) + 
  geom_point(size = 1) + 
  scale_colour_gradientn(colours = c("black", "purple", "red", "orange", "yellow", "lightgreen"), 
                         name = expression(italic(r)), 
                         breaks = pretty_breaks(n = 3)) +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  labs(x = "Nº of model runs", y = "$k$") + 
  facet_wrap(~variable, ncol = 1, labeller = labeller(variable = supp.labs)) +
  theme_AP() + 
  theme(legend.position = "none", 
        strip.background = element_rect(fill = "white"))

b <- melt(final.output, measure.vars = discrepancy_methods) %>%
  ggplot(., aes(reorder(variable, -value), value)) +
  geom_boxplot() + 
  coord_flip() +
  labs(x = "", y = "$r$") + 
  scale_x_discrete(position = "top", 
                   labels = supp.labs) +
  theme_AP()

b

legend <- get_legend(a + theme(legend.position = "top", 
                               legend.margin=margin(0, 0, 0, 0),
                               legend.box.margin=margin(-7,-7,-7,-7)))
bottom <- plot_grid(a, b, ncol = 2, labels = "auto")

plot_grid(legend, bottom, ncol = 1, rel_heights = c(0.08, 0.9))


# SOME STATISTICS -------------------------------------------------------------

final.stat <- melt(final.output, measure.vars = discrepancy_methods) 

final.stat[, .(mean = mean(value), median = median(value), sd = sd(value), 
        max = max(value), min = min(value)), variable]

final.stat[, .(quantile = quantile(value)), variable]







## ----session_information---------------------------------------------------------

# SESSION INFORMATION --------------------------------------------------------------

sessionInfo()

## Return the machine CPU
cat("Machine:     "); print(get_cpu()$model_name)

## Return number of true cores
cat("Num cores:   "); print(detectCores(logical = FALSE))

## Return number of threads
cat("Num threads: "); print(detectCores(logical = TRUE))

## Return the machine RAM
cat("RAM:         "); print (get_ram()); cat("\n")










prove <- final.mat[1, ]
k <- prove[["k"]]
params <- paste("X", 1:k, sep = "")
tau <- prove[["tau"]]
epsilon <- prove[["epsilon"]]
phi <- prove[["phi"]]
cost.discrepancy <- prove[["cost.discrepancy"]]
base.sample.size <- prove[["base.sample.size"]]




if (tau == 1) {
  
  type <- "R"
  
} else if (tau == 2) {
  
  type <- "QRN"
}
set.seed(epsilon)
discrepancy.mat <- sobol_matrices(matrices = "A", N = cost.discrepancy, 
                                  params = params, 
                                  type = type)

set.seed(epsilon)
jansen.mat <- sobol_matrices(matrices = c("A", "AB"), N = base.sample.size, 
                             params = params, 
                             type = type)


all.matrices <- rbind(discrepancy.mat, jansen.mat)
set.seed(epsilon)
transformed.all.matrices <- random_distributions(X = all.matrices, phi = phi)

y <- sensobol::metafunction(data = transformed.all.matrices, epsilon = epsilon)

type <- c("symmetric", "star", "L2", "centered")

discrepancy.value <- lapply(type, function(x) 
  discrepancy(mat = all.matrices[1:cost.discrepancy, ], 
              y = y[1:cost.discrepancy], params = params, type = x))


jansen.value <- jansen_ti(d = y[(cost.discrepancy + 1):length(y)],
                          N = base.sample.size, params = params)

savage.discrepancy <- lapply(discrepancy.value, savage_scores)
jansen.discrepancy <- savage_scores(jansen.value)

out <- lapply(savage.discrepancy, function(x) 
  cor(x, jansen.discrepancy))

plot_scatter(data = all.matrices[1:cost.discrepancy, ], 
             N = cost.discrepancy, Y = y[1:cost.discrepancy], params = params)






da <- data.table("jansen" = jansen.value, 
           "discrepancy" = discrepancy.value[[1]]) %>%
  .[, rank.jansen:= rank(-jansen)] %>%
  .[, rank.discrepancy:= rank(-discrepancy)] %>%
  .[, parameter:= paste("X", 1:nrow(.), sep = "")]


da




mat




params <- paste("X", 1:k, sep = "")

if (tau == 1) {
  
  type <- "R"
  
} else if (tau == 2) {
  
  type <- "QRN"
}
set.seed(epsilon)
discrepancy.mat <- sobol_matrices(matrices = "A", N = cost.discrepancy, 
                                  params = params, 
                                  type = type)

set.seed(epsilon)
jansen.mat <- sobol_matrices(matrices = c("A", "AB"), N = base.sample.size, 
                             params = params, 
                             type = type)

set.seed(epsilon)
all.matrices <- random_distributions(X = rbind(discrepancy.mat, jansen.mat), phi = phi)

# all.matrices <- rbind(discrepancy.mat, jansen.mat)

y <- sensobol::metafunction(data = all.matrices, epsilon = epsilon)

type <- c("symmetric", "star", "L2", "centered")

discrepancy.value <- lapply(type, function(x) 
  discrepancy(mat = all.matrices[1:cost.discrepancy, ], 
              y = y[1:cost.discrepancy], params = params, type = x))


jansen.value <- jansen_ti(d = y[(cost.discrepancy + 1):length(y)],
                          N = base.sample.size, params = params)

savage.discrepancy <- lapply(discrepancy.value, savage_scores)
jansen.discrepancy <- savage_scores(jansen.value)

out <- lapply(savage.discrepancy, function(x) 
  cor(x, jansen.discrepancy))
