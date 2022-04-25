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
               "cowplot", "pcaPP"))


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
                   "L2wraparound_functions.cpp", "L2modified_functions.cpp")

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
    
  } else if (type == "modified") {
    
    P <- 3 - X^2
    s1 <- DisM2_Rowprod(t(P), dimension)
    s2 <- DisM2_Crossprod(c(t(X)), dimension)
    R <- sqrt(((4/3)^dimension) - (((2^(1 - dimension))/n) * s1) + ((1/n^2) * s2))
    
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
fun_list <- list(bratley1992_Fun, oakley_Fun, sensitivity::morris.fun)
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
discrepancy_methods <- c("symmetric", "star", "L2", "centered", 
                         "wraparound", "modified")

check_discrepancy <- function(type) {
  
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
      output[[i]][[j]] <- lapply(sample.sizes, function(x) 
        cor_fun(N = x, params = params, replicas = replicas, model_fun = fun_list[[i]], 
                true_scores = savage.rank[[i]], approach = j, type = type,
                scrambling = 1))
      
      names(output[[i]][[j]]) <- final.sample.sizes
      
    }
  }
  return(output)
}

output <- mclapply(discrepancy_methods, function(type) 
  check_discrepancy(type = type), mc.cores = detectCores() * 0.75)

names(output) <- discrepancy_methods

# Arrange output -------------------
results <- lapply(output, function(x) lapply(x, function(y) 
  lapply(y, function(z) lapply(z, function(w) data.table(w))))) %>%
  lapply(., function(x) lapply(x, function(y) lapply(y, function(z)  
    rbindlist(z, idcol = "Model.runs")))) %>%
  lapply(., function(x) lapply(x, function(y) rbindlist(y, idcol = "Approach"))) %>%
  lapply(., function(x) rbindlist(x, idcol = "Function")) %>%
  rbindlist(., idcol = "Method") %>%
  .[, Model.runs:= as.numeric(Model.runs)] %>%
  .[, Approach:= ifelse(Approach == "discrepancy", "Discrepancy", "Jansen estimator")] %>%
  na.omit()

## ----plot, dependson="model", warning=FALSE, fig.height=4, fig.cap="Correlation between the Savage scores obtained with $T_i$ (Jansen estimator) and those obtained with discrepancy."----

# Plot the results ------------------------------------------------------------

# New facet label names for supp variable
supp.labs <- c("Symmetric $L2$ discrepancy", "Jansen estimator")
names(supp.labs) <- approach.vec

a <- ggplot(results[Approach == "Jansen estimator"], aes(Model.runs, w, group = Model.runs)) +
  geom_boxplot(position = "identity", outlier.size = 0.5) +
  labs(x = "Nº of model runs", y = "Correlation") +
  facet_grid(~Function, 
             scale = "free_x", space = "free") +
  theme_AP()

b <- ggplot(results[Approach == "Discrepancy"], aes(Model.runs, w, group = Model.runs)) +
  geom_boxplot(position = "identity", outlier.size = 0.5) +
  labs(x = "Nº of model runs", y = "Correlation") +
  facet_grid(Method~Function, 
             scale = "free_x", space = "free") +
  theme_AP()

dev.off()
da <- list(a, b)
plot_grid(plotlist = da, ncol = 1, labs = c("a", "b"))
ggpubr::ggarrange(plotlist = list(a, b), ncol = 1, labs = "auto")
plot_grid(a, b, ncol = 1, labs = "auto", rel_heights = 0.5, 1)

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

ggplot(data.frame(x = runif(100)), aes(x)) +
  map(1:length(function_list), function(nn) {
    stat_function(fun = function_list[[nn]], 
                  geom = "line", 
                  aes_(color = factor(names(function_list[nn]))))
  }) + 
  labs(color= "Function", linetype = "Function", 
       x = expression(italic(x)), 
       y = expression(italic(y))) +
  scale_color_discrete(labels = c("$f_1(x) = x^3$", 
                                  "$f_2(x) = 1 \\hspace{1mm} \\mbox{if} \\hspace{1mm} x > \\frac{1}{2}, \\hspace{1mm} \\mbox{otherwise} \\hspace{1mm} 0$", 
                                  "$f_3(x) = \\frac{(e^x - 1)}{e-1}$", 
                                  "$f_4(x) = (10-\\frac{1}{1.1})^{-1}(x + 0.1)^{-1}$", 
                                  "$f_5(x) = x$", 
                                  "$f_6(x) = 0$", 
                                  "$f_7(x) = 4(x - 0.5)^2$", 
                                  "$f_8(x) = 2 - 0.2 \\cos(7 \\pi x)$", 
                                  "$f_9(x) = \\frac{\\sin(2 \\pi x)}{2}$", 
                                  "$f_{10}(x) = (-1)^ {|4x|} [0.125- \\mbox{mod}(x, 0.25)] + 0.125$", 
                                  "$f_{11}(x) = (-1)^ {|32x|} [0.0325-\\mbox{mod}(x, 0.0325)] + 0.0325$", 
                                  "$f_{12}(x) = x^2$", 
                                  "$f_{13}(x) = \\cos(x)$")) +
  theme_AP() + 
  theme(legend.text.align = 0)


k_epsilon <- list(c(3, 2), c(12, 6))
out <- list()

for(i in 1:length(k_epsilon)) {
  params <- paste("$x_", 1:k_epsilon[[i]][[1]], "$", sep = "")
  mat <- sobol_matrices(N = N, params = params)
  y <- metafunction(mat, epsilon = k_epsilon[[i]][[2]])
  out[[i]] <- plot_scatter(data = mat, N = N, Y = y, params = params) + 
    scale_x_continuous(breaks = pretty_breaks(n = 3)) +
    theme(plot.margin = unit(c(0, 0.1, 0, 0), "cm")) + 
    labs(x = "", y = "$y$")

}


N <- 2^9
params <- paste("X", 1:13, sep = "")
mat <- sobol_matrices(N = N, params = params)
y <- 

function_names <- names(function_list)

# Compute model output first order effects
mat.y <- sapply(seq_along(function_names), function(x)
  function_list[[function_names[x]]](mat[, x]))

# Compute first-order effects
y1 <- Rfast::rowsums(mat.y)

plot_scatter(data = mat, N = N, Y = y1, params = params)



scatters_plot <- plot_grid(out[[1]] + facet_wrap(~variable, ncol = 1), 
          out[[2]] + facet_wrap(~variable, ncol = 4) + labs(y = ""), 
          ncol = 2, labels = "auto", rel_widths = c(0.24, 0.76))

scatters_plot







N <- 2^9
k <- 12
params <- paste("X", 1:k, sep = "")
mat <- sobol_matrices(N = N, params = params)
y <- metafunction(mat, epsilon = 3)
plot_scatter(data = mat, N = N, Y = y, params = params) + 
  scale_x_continuous(breaks = pretty_breaks(n = 3))

## ----model_fun, dependson=c("metamodel_fun", "savage_fun", "discrepancy_fun", "jansen_fun")----

# DEFINE MODEL ----------------------------------------------------------------

model_fun <- function(tau, epsilon, base.sample.size, 
                      cost.discrepancy, phi, k) {
  
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
  
  type <- c("symmetric", "star", "L2", "centered", "wraparound", "modified")
  
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

N <- 2^4
params <- c("epsilon", "phi", "k", "tau",  "base.sample.size")
mat <- sobol_matrices(matrices = "A", N = N, params = params)

# Define distributions --------------------------

mat[, "epsilon"] <- floor(qunif(mat[, "epsilon"], 1, 200))
mat[, "phi"] <- floor(mat[, "phi"] * 8) + 1
mat[, "k"] <- floor(qunif(mat[, "k"], 3, 50))
mat[, "tau"] <- floor(mat[, "tau"] * 2) + 1
mat[, "base.sample.size"] <- floor(qunif(mat[, "base.sample.size"], 10, 50))

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
discrepancy_methods <- c("symmetric", "star", "L2", "centered", 
                         "wraparound", "modified")
colnames(output) <- discrepancy_methods

final.output <- data.table(cbind(final.mat, output))

# Compute mean and median
melt(final.output, measure.vars = discrepancy_methods) %>%
  .[, .(mean = mean(value), median = median(value)), variable]

## ----plot_uncertainty, dependson="run_model_fun", warning = FALSE, fig.height=3, fig.width=3----

# PLOT UNCERTAINTY ------------------------------------------------------------

# New facet label names for supp variable
supp.labs <- c("Symmetric", "Star", "$L_2$", "Centered", "Wrap-around", "Modified")
names(supp.labs) <- discrepancy_methods

a <- melt(final.output, measure.vars = discrepancy_methods) %>%
  .[, variable:= factor(variable, levels =   c("symmetric", "wraparound", 
                                               "centered", "L2",
                                               "star","modified"))] %>%
  ggplot(., aes(cost.discrepancy, k, color = value)) + 
  geom_point(size = 1) + 
  scale_colour_gradientn(colours = c("black", "purple", "red", "orange", "yellow", "lightgreen"), 
                         name = expression(italic(r)), 
                         breaks = pretty_breaks(n = 3)) +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  labs(x = "Nº of model runs", y = "$k$") + 
  facet_wrap(~variable, ncol = 6, labeller = labeller(variable = supp.labs)) +
  theme_AP() + 
  theme(legend.position = "none", 
        strip.background = element_rect(fill = "white"))

b <- melt(final.output, measure.vars = discrepancy_methods) %>%
  ggplot(., aes(reorder(variable, -value), value)) +
  geom_boxplot() + 
  labs(x = "", y = "$r$") + 
  scale_x_discrete(position = "bottom", 
                   labels = supp.labs) +
  theme_AP()

b

legend <- get_legend(a + theme(legend.position = "top", 
                               legend.margin=margin(0, 0, 0, 0),
                               legend.box.margin=margin(-7,-7,-7,-7)))
bottom <- plot_grid(a, b, ncol = 1, labels = c("c", "d"))

da <- plot_grid(legend, bottom, ncol = 1, rel_heights = c(0.15, 0.85))
plot_grid(scatters_plot, da, ncol = 1, rel_heights = c(0.45, 0.55))



final.output <- fread("final.output.csv")

final.output[, ratio:= cost.discrepancy / k] %>%
  melt(., measure.vars = discrepancy_methods) %>%
  ggplot(., aes(ratio, value)) +
  geom_point(alpha = 0.1, size = 0.2) +
  facet_wrap(~variable, 
             ncol = 1) +
  geom_smooth() +
  labs(x = expression(italic(N[t]/k)), 
       y = expression(italic(r))) +
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  theme_AP() +
  theme(strip.background = element_rect(fill = "white"))
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


final.stat <- melt(final.output, measure.vars = discrepancy_methods) 

final.stat[, .(mean = mean(value), median = median(value), sd = sd(value), 
               max = max(value), min = min(value)), variable] %>%
  .[order(-median)]


###################################

# NEW FUNCTIONS ------------------------------------

f1_fun <- function(x) 10 * x[, 1] + 0.2 * x[, 2]^3

f2_fun <- function(x) 2 * x[, 1] - x[, 2]^2

f3_fun <- function(x) x[, 1]^2 + x[, 2]^4 + x[, 1] * x[, 2] + x[, 2] * x[, 3]^4

f4_fun <- function(x) 0.2 * exp(x[, 1] - 3) + 2.2 * abs(x[, 2]) + 1.3 * x[, 2]^6 -
  2 * x[, 2]^2 - 0.5 * x[, 2]^4 - 0.5 * x[, 1]^4 + 2.5 * x[, 1]^2 + 0.7 * x[, 1]^3 +
  3 / ((8 * x[, 1] - 2)^2 + (5 * x[, 2] - 3)^2 + 1) + sin(5 * x[, 1]) * cos(3 * x[, 1]^2)

f6_fun <- function(x) 0.2 * exp(x[, 1] + 2 * x[, 2])

fun_vec <- paste("f", 1:4, "_fun", sep = "")

f_list <- list(f1_fun, f2_fun, f3_fun, f4_fun)
names(f_list) <- fun_vec

# RUN FUNCTIONS ---------------------------------------------------------------

output <- ind <- list()

for (i in names(f_list)) {
  
  if(i == "f3_fun") {
    
    k <- 3
    
  } else {
    
    k <- 2
  }
  params <- paste("$x_", 1:k, "$", sep = "")
  N <- 2^9
  mat <- sobol_matrices(N = N, params = params, scrambling = 1)
  y <- f_list[[i]](mat)
  y <- rescale_fun(y)
  ind[[i]] <- sobol_indices(Y = y, N = N, params = params)
  output[[i]] <- plot_scatter(data = mat, N = N, Y = y, params = params) + 
    labs(x = "$x$", y = "$y$") +
    scale_x_continuous(breaks = pretty_breaks(n = 3)) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
}

plot_list <- list(output[[1]] + labs(x = "", y = "$y$"), 
                  output[[2]] + labs(x = "", y = "$y$"), 
                  output[[4]])

plot_grid(plotlist = plot_list, ncol = 1, labels = "auto")


plot_grid(bottom, output[[3]] + 
            facet_wrap(~variable, ncol = 1, scales = "free_x"), ncol = 2, labels = c("", "c"))


##########################3

# PLOT TOTAL-ORDER INDEX -----------------------------------------------------

plot_ti <- function(ind, params) {
  
  out <- data.table(ind) %>%
    .[, parameters:= params] %>%
    ggplot(., aes(parameters, ind)) +
    geom_bar(stat = "identity", fill = "#00BFC4", color = "black") + 
    labs(x = "", y = "$T_i$") +
    theme_AP()
  
  return(out)
}

# FUNCTION TO DISCRETIZE THE BRATLEY ET AL. FUNCTION --------------------------

discretization_fun <- function(N, k, epsilon, plot.scatterplot = FALSE, 
                               plot.jansen = FALSE, savage.scores = TRUE) {

  params <- paste("$x_", 1:k, "$", sep = "")
  mat <- sobol_matrices(N = N, params = params, matrices = c("A", "AB"))
  mat2 <- mat # mat2 = matrix for discretization
  
  # Determine how many parameters to discretize (n.discrete), 
  # which parameters to discretize (n.parameters), and the 
  # number of levels (n.levels)
  set.seed(epsilon)
  n.discrete <- sample(1:k, size = 1)
  set.seed(epsilon)
  n.parameters <- sample(1:k, size = n.discrete)
  set.seed(epsilon)
  n.levels <- sample(1:10, size = 1)
  
  # Discretization of the sampled parameters
  mat2[, n.parameters] <- floor(mat2[, n.parameters] * n.levels) + 1
  
  # Compute output and sobol indices (Jansen estimator)
  y <- bratley1992_Fun(mat2)
  ind2 <- jansen_ti(d = y, N = N, params = params)
  
  # Conditional
  if (plot.scatterplot == TRUE) {
    out <- plot_scatter(data = mat2, N = N, Y = y, params = params)
  } 
  if (plot.jansen == TRUE) {
    out <- plot_ti(ind = ind2, params = params)
  }
  if(savage.scores == TRUE) {
    sobol_ranks <- savage_scores(ind2)
    
    # Compute discrepancy
    type <- c("symmetric", "star", "L2", "centered", "wraparound", "modified")
    discrepancy.out <- lapply(type, function(x) 
      discrepancy(mat = mat[1:N, ], y = y[1:N], params = params, type = x))
    
    discrepancy.ranks <- lapply(discrepancy.out, savage_scores)
    # Correlation on savage scores
    out <- lapply(discrepancy.ranks, function(x) cor(sobol_ranks, x))
  }
  return(out)
}

# CREATE SAMPLE MATRIX -------------------------------------------------------

N <- 2^5
params <- c("N", "k", "epsilon")
mat <- sobol_matrices(matrices = "A", N = N, params = params)

# Define distributions --------------------------

mat[, "epsilon"] <- floor(qunif(mat[, "epsilon"], 1, N))
mat[, "N"] <- floor(qunif(mat[, "N"], 1000, 2000))
mat[, "k"] <- floor(qunif(mat[, "k"], 6, 6))

# RUN MODEL --------------------------------------------------------------------

y.discrete <- mclapply(1:nrow(mat), function(i) {
  discretization_fun(N = mat[i, "N"],
                     k = mat[i, "k"],
                     epsilon = mat[i, "epsilon"], 
                     plot.scatterplot = FALSE, 
                     plot.jansen = FALSE, 
                     savage.scores = TRUE)}, 
  mc.cores = floor(detectCores() * 0.75))

# RUN MODEL TO GET SCATTERPLOTS -----------------------------------------------

y.scatters <- mclapply(1:nrow(mat[1:10, ]), function(i) {
  discretization_fun(N = mat[i, "N"],
                     k = mat[i, "k"],
                     epsilon = mat[i, "epsilon"], 
                     plot.scatterplot = TRUE, 
                     plot.jansen = FALSE, 
                     savage.scores = FALSE)}, 
  mc.cores = floor(detectCores() * 0.75))

# RUN MODEL TO GET JANSEN TI --------------------------------------------------

y.indices <- mclapply(1:nrow(mat[1:10, ]), function(i) {
  discretization_fun(N = mat[i, "N"],
                     k = mat[i, "k"],
                     epsilon = mat[i, "epsilon"], 
                     plot.scatterplot = FALSE, 
                     plot.jansen = TRUE, 
                     savage.scores = FALSE)}, 
  mc.cores = floor(detectCores() * 0.75))

# ARRANGE OUTPUT ---------------------------------------------------------------

y.discrete <- lapply(y.discrete, unlist)
output <- data.table(do.call(rbind, y.discrete))
discrepancy_methods <- c("symmetric", "star", "L2", "centered", 
                         "wraparound", "modified")
colnames(output) <- discrepancy_methods

final.output.discrete <- data.table(cbind(mat, output))

melt(final.output.discrete, measure.vars = discrepancy_methods) %>%
  .[, .(mean = mean(value, na.rm = TRUE), 
        median = median(value, na.rm = TRUE)), variable] %>%
  .[order(-median)]

# PLOT UNCERTAINTY -------------------------------------------------------------

boxplots.discrete <- melt(final.output.discrete, measure.vars = discrepancy_methods) %>%
  ggplot(., aes(reorder(variable, -value), value)) +
  geom_boxplot() + 
  scale_x_discrete(guide = guide_axis(n.dodge = 2), 
                   labels = supp.labs) +
  labs(x = "", y = "$r$") + 
  theme_AP() +
  theme(plot.margin = unit(c(0, 0.1, 0, 0), "cm"))

boxplots.discrete

# PLOT ORIGINAL BRATLEY ET AL. FUNCTION ---------------------------------------

# Compute Ti -----------------------------
N <- 2^9
k <- 13
params <- paste("$x_", 1:k,"$", sep = "")
matrices <- c("A", "AB")
mat <- sobol_matrices(N = N, params = params, matrices = matrices)
y <- bratley1992_Fun(mat)
ind <- jansen_ti(d = y, N = N, params = params)


bratley_scatter <- plot_scatter(data = mat, N = N, Y = y, params = params) +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) + 
  scale_y_continuous(breaks = pretty_breaks(n = 3)) + 
  facet_wrap(~variable, ncol = 6) +
  labs(x = "$x$", y = "$y$") +
  theme(plot.margin = unit(c(0, 0.1, 0, 0), "cm"))
  
bratley_ti <- plot_ti(ind = ind, params = params) + 
  theme(plot.margin = unit(c(0, 0.1, 0, 0), "cm")) + 
  labs(x = "", y = "$T_i$") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

bratley_plots <- plot_grid(bratley_scatter, bratley_ti, labels = c("a", "b"), 
                           ncol = 2, rel_widths = c(0.7, 0.3))
bratley_plots

plot_grid(boxplots.discrete, bratley_plots, labels = c("a", ""), ncol = 1, 
          rel_heights = c(0.4, 0.6))



selected.plots <- 2

plot_grid(y.scatters[[selected.plots]], y.indices[[selected.plots]], 
          ncol = 1)



#########################################################

N <- 2^11
k <- 6
params <- paste("X", 1:k, sep = "")
R <- 10^2
matrices <- c("A", "B", "AB")

mat <- sobol_matrices(N = N, params = params)
y <- bratley1992_Fun(mat)
ind <- sobol_indices(Y = y, N = N, params = params, matrices = matrices, 
                     R = R, boot = TRUE)
plot(ind)

sub.sample <- seq(500, N, 50)
sobol_convergence(matrices = matrices, Y = y, N = N, sub.sample = sub.sample, 
                  params = params, first = "saltelli", total = "jansen", order = "first", 
                  plot.order = "first")








# ADJUSTED FOR TOTAL NUMBER OF MODEL RUNS
#########################################

cost.discrepancy <- N * (k + 1)
params <- paste("X", 1:k, sep = "")
mat.jansen <- sobol_matrices(N = N, params = params, matrices = c("A", "AB"))
# matrix for discrepancy
mat.discrepancy <- sobol_matrices(N = cost.discrepancy, params = params, matrices = "A")
mat.discrepancy2 <- mat.discrepancy

# Determine how many parameters to discretize (n.discrete), 
# which parameters to discretize (n.parameters), and the 
# number of levels (n.levels)
set.seed(epsilon)
n.discrete <- sample(1:k, size = 1)
set.seed(epsilon)
n.parameters <- sample(1:k, size = n.discrete)
set.seed(epsilon)
n.levels <- sample(1:10, size = 1)

# Discretization of the sampled parameters
mat.jansen[, n.parameters] <- floor(mat.jansen[, n.parameters] * n.levels) + 1
mat.discrepancy[, n.parameters] <- floor(mat.discrepancy[, n.parameters] * n.levels) + 1

all.matrices <- rbind(mat.discrepancy, mat.jansen)

# Compute output and sobol indices (Jansen estimator)
y <- bratley1992_Fun(all.matrices)
ind <- jansen_ti(d = y[(cost.discrepancy + 1):length(y)], N = N, params = params)
sobol_ranks <- savage_scores(ind)

# Compute discrepancy
type <- c("symmetric", "star", "L2", "centered", "wraparound", "modified")
discrepancy.out <- lapply(type, function(x) 
  discrepancy(mat = mat.discrepancy2, y = y[1:cost.discrepancy], params = params, type = x))

discrepancy.ranks <- lapply(discrepancy.out, savage_scores)

# Correlation on savage scores
out <- lapply(discrepancy.ranks, function(x) cor(sobol_ranks, x))







k <- 6
N <- 2^9
epsilon <- 2

cost.discrepancy <- N * (k + 1)
params <- paste("X", 1:k, sep = "")
mat.jansen <- sobol_matrices(N = N, params = params, matrices = c("A", "AB"))
# matrix for discrepancy
mat.discrepancy <- sobol_matrices(N = cost.discrepancy, params = params, matrices = "A")

# Determine how many parameters to discretize (n.discrete), 
# which parameters to discretize (n.parameters), and the 
# number of levels (n.levels)
set.seed(epsilon)
n.discrete <- sample(1:k, size = 1)
set.seed(epsilon)
n.parameters <- sample(1:k, size = n.discrete)
set.seed(epsilon)
n.levels <- sample(1:10, size = 1)

# Discretization of the sampled parameters
mat.jansen[, n.parameters] <- floor(mat.jansen[, n.parameters] * n.levels) + 1
mat.discrepancy[, n.parameters] <- floor(mat.discrepancy[, n.parameters] * n.levels) + 1

# Compute output and sobol indices (Jansen estimator)
y.jansen <- bratley1992_Fun(mat.jansen)
y.discrepancy <- bratley1992_Fun(mat.discrepancy)
ind2 <- jansen_ti(d = y.jansen, N = N, params = params)


ind2 <- sobol_indices(Y = y.jansen, N = N, params = params)
sobol_ranks <- savage_scores(ind2$results[sensitivity == "Ti", original])

# Compute discrepancy
type <- c("symmetric", "star", "L2", "centered", "wraparound", "modified")
discrepancy.out <- lapply(type, function(x) 
  discrepancy(mat = mat.discrepancy, y = y.discrepancy, params = params, type = x))

discrepancy.ranks <- lapply(discrepancy.out, savage_scores)

# Correlation on savage scores
out <- lapply(discrepancy.ranks, function(x) cor(sobol_ranks, x))

















































# Define settings (test with k = 10)
N <- 2^11
params <- paste("X", 1:6, sep = "")

# Create sample matrix
mat <- sobol_matrices(N = N, params = params)

# Compute Bratley et al. (1992) function
y <- bratley1992_Fun(mat)
ind <- sobol_indices(Y = y, N = N, params = params)
plot(ind)
plot_scatter(data = mat, N = N, Y = y, params = params)













all_funs <- function(matrices, N, params, name_fun) {
  mat <- sobol_matrices(matrices = matrices, N = N, params = params)
  y <- f_list[[name_fun]](mat)
  ind <- sobol_indices(matrices = matrices, Y = y, N = N, params = params)
  p <- plot_scatter(data = mat, N = N, Y = y, params = params)
  out <- list(ind, p)
  return(out)
}



N <- 2^10
k <- 2
params <- paste("X", 1:k, sep = "")
name_fun <- "f2_fun"
matrices <- c("A", "B", "AB")

out <- all_funs(matrices = matrices, N = N, params = params, name_fun = name_fun)
out












for (i in names(f_list)) {
  if (i == "f3_fun") {
    k <- 3
  } else {
    k <- 2
  }  
  
}





sobol_Fun()

N <- 2^10
k <- 2
params <- paste("X", 1:k, sep = "")
mat <- sobol_matrices(N = N, params = params, scrambling = 1)

y <- f4_fun(mat)
ind <- sobol_indices(Y = y, N = N, params = params)
ind
plot(ind)

plot_scatter(data = mat, N = N, Y = y, params = params)

discrepancy(mat = mat, y = y, params = params, type = "symmetric")

#################



kriging_fun <- function(x) 1 + exp(-2 * ((x[, 1] - 1)^2 + x[, 2]^2) - 0.5 * 
  (x[, 3]^2 + x[, 4]^2)) + exp(-2 * (x[, 1]^2 + (x[, 2] -1)^2) - 0.5 * 
                                 (x[, 3]^2 + x[, 4]^2))


N <- 2^10
k <- 4
params <- paste("X", 1:k, sep = "")
mat <- sobol_matrices(N = N, params = params)
y <- kriging_fun(mat)
ind <- sobol_indices(Y = y, N = N, params = params)
plot(ind)


