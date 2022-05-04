## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, dev = "tikz", cache = TRUE)


## ---- results="hide", message=FALSE, warning=FALSE, cache=FALSE----------

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
               "cowplot", "wesanderson"))


## ----source_cpp----------------------------------------------------------

# CPP CODE ---------------------------------------------------------------------

# Source cpp code -------------------------------
cpp_functions <- c("cpp_functions.cpp", "L2star_functions.cpp", 
                   "L2_functions.cpp", "L2centered_functions.cpp", 
                   "L2wraparound_functions.cpp", "L2modified_functions.cpp")

for(i in 1:length(cpp_functions)) {
  Rcpp::sourceCpp(cpp_functions[i])
}


## ----savage_fun----------------------------------------------------------

# SAVAGE SCORES FUNCTION -------------------------------------------------------

savage_scores <- function(x) {
  true.ranks <- rank(-x)
  p <- sort(1 / true.ranks)
  mat <- matrix(rep(p, length(p)), nrow = length(p), byrow = TRUE)
  mat[upper.tri(mat)] <- 0
  out <- sort(rowSums(mat), decreasing = TRUE)[true.ranks]
  return(out)
}


## ----discrepancy_fun, dependson="savage_scores"--------------------------

# DISCREPANCY FUNCTION --------------------------------------------------------

# This is based on the function discrepancyCriteria_cpp of
# the sensitivity package

# Function to rescale ---------------------------
rescale_fun <- function(x) (x - min(x)) / (max(x) - min(x))

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


## ----jansen_fun----------------------------------------------------------

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


## ----fun_list------------------------------------------------------------

# LIST OF FUNCTIONS -----------------------------------------------------------

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
  N <- 2^7
  mat <- sobol_matrices(N = N, params = params, scrambling = 1)
  y <- f_list[[i]](mat)
  y <- rescale_fun(y)
  ind[[i]] <- sobol_indices(Y = y, N = N, params = params)
  output[[i]] <- plot_scatter(data = mat, N = N, Y = y, params = params) + 
    labs(x = "$x$", y = "$y$") +
    scale_x_continuous(breaks = pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = pretty_breaks(n = 3)) +
    theme(plot.margin = unit(c(0, 0.2, 0, 0), "cm"))
}

ind


## ----plot_fun_list, dependson="fun_list", fig.height=4, fig.width=2.5----

# PLOT LIST OF FUNCTIONS -------------------------------

plot_list <- list(output[[1]] + labs(x = "", y = "$y$"), 
                  output[[2]] + labs(x = "", y = "$y$"), 
                  output[[4]])

scat_plot <- plot_grid(plotlist = plot_list, ncol = 1, labels = "auto")
scat_plot


## ----metamodel_fun-------------------------------------------------------

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


## ----plot_metafunction, dependson="metamodel_fun", fig.width=5.5, fig.height=3.3----

# PLOT METAFUNCTION -----------------------------------------------------------

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


## ----scatters_plot, fig.height=3, fig.width=5.5--------------------------

# PLOT SCATTERS FROM METAFUNCTION ---------------------------------------------

k_epsilon <- list(c(3, 2), c(12, 6))
out <- list()

for(i in 1:length(k_epsilon)) {
  N <- 2^9
  params <- paste("$x_{", 1:k_epsilon[[i]][[1]], "}$", sep = "")
  mat <- sobol_matrices(N = N, params = params)
  y <- metafunction(mat, epsilon = k_epsilon[[i]][[2]])
  out[[i]] <- plot_scatter(data = mat, N = N, Y = y, params = params) + 
    scale_x_continuous(breaks = pretty_breaks(n = 3)) +
    theme(plot.margin = unit(c(0, 0.1, 0, 0), "cm")) + 
    labs(x = "", y = "$y$")

}

scatters_plot <- plot_grid(out[[1]] + facet_wrap(~variable, ncol = 1), 
          out[[2]] + facet_wrap(~variable, ncol = 4) + labs(y = ""), 
          ncol = 2, labels = "auto", rel_widths = c(0.24, 0.76))

scatters_plot


## ----model_fun, dependson=c("metamodel_fun", "savage_fun", "discrepancy_fun", "jansen_fun")----

# DEFINE MODEL -----------------------------------------------------------------

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


## ----matrix--------------------------------------------------------------

# CREATE SAMPLE MATRIX --------------------------------------------------------

N <- 2^10
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


## ----run_model_fun, dependson=c("model_fun", "matrix")-------------------

# RUN MODEL  -----------------------------------------------------------------

y <- mclapply(1:nrow(final.mat), function(i) {
  model_fun(tau = final.mat[i, "tau"],
            epsilon = final.mat[i, "epsilon"], 
            base.sample.size = final.mat[i, "base.sample.size"], 
            cost.discrepancy = final.mat[i, "cost.discrepancy"], 
            phi = final.mat[i, "phi"], 
            k = final.mat[i, "k"])}, 
  mc.cores = floor(detectCores() * 0.75))


## ----arrange_output, dependson=c("run_model_fun", "matrix")--------------

# ARRANGE DATA ----------------------------------------------------------------

y <- lapply(y, unlist)
output <- data.table(do.call(rbind, y))
discrepancy_methods <- c("symmetric", "star", "L2", "centered", 
                         "wraparound", "modified")
colnames(output) <- discrepancy_methods

final.output <- data.table(cbind(final.mat, output))

# Write model output
fwrite(final.output, "final.output.csv")


## ----discretization_fun, dependson = c("jansen_fun", "discrepancy_fun", "savage_fun")----

# DISCRETIZATION FUNCTIONS -----------------------------------------------------

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
  n.discrete <- sample(1:ceiling(k / 2), size = 1)
  set.seed(epsilon)
  n.parameters <- sample(1:k, size = n.discrete)
  set.seed(epsilon)
  n.levels <- sample(2:5, size = 1)
  
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


## ----discrete_sims, dependson=c("jansen_fun", "discrepancy_fun", "savage_fun", "discretization_fun")----

# CREATE SAMPLE MATRIX -------------------------------------------------------

N <- 2^7
params <- c("N", "k", "epsilon")
mat <- sobol_matrices(matrices = "A", N = N, params = params)

# Define distributions --------------------------

mat[, "epsilon"] <- floor(qunif(mat[, "epsilon"], 1, N))
mat[, "N"] <- floor(qunif(mat[, "N"], 5000, 5000))
mat[, "k"] <- floor(qunif(mat[, "k"], 6, 50))

# RUN MODEL --------------------------------------------------------------------

y.discrete <- mclapply(1:nrow(mat), function(i) {
  discretization_fun(N = mat[i, "N"],
                     k = mat[i, "k"],
                     epsilon = mat[i, "epsilon"], 
                     plot.scatterplot = FALSE, 
                     plot.jansen = FALSE, 
                     savage.scores = TRUE)}, 
  mc.cores = floor(detectCores() * 0.75))


## ----arrange_discrete, dependson="discrete_sims"-------------------------

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


## ----discretization_scatter, dependson="discrete_sim"--------------------

# RUN MODEL TO GET SCATTERPLOTS -----------------------------------------------

N.scatter <- 2^8
N.indices <- 2^10

y.scatters <- mclapply(1:nrow(mat[1:10, ]), function(i) {
  discretization_fun(N = N.scatter,
                     k = mat[i, "k"],
                     epsilon = mat[i, "epsilon"], 
                     plot.scatterplot = TRUE, 
                     plot.jansen = FALSE, 
                     savage.scores = FALSE)}, 
  mc.cores = floor(detectCores() * 0.75))

# RUN MODEL TO GET JANSEN TI --------------------------------------------------

y.indices <- mclapply(1:nrow(mat[1:10, ]), function(i) {
  discretization_fun(N = N.indices,
                     k = mat[i, "k"],
                     epsilon = mat[i, "epsilon"], 
                     plot.scatterplot = FALSE, 
                     plot.jansen = TRUE, 
                     savage.scores = FALSE)}, 
  mc.cores = floor(detectCores() * 0.75))


## ----bratley_fun, dependson="discrepancy_fun", fig.height=1.3, fig.width=6.4----

# PLOT ORIGINAL BRATLEY ET AL. FUNCTION ---------------------------------------

# Settings -----------------------------
N <- 2^7
k <- 6
params <- paste("$x_", 1:k,"$", sep = "")
matrices <- c("A", "AB")

# Compute -----------------------------
mat <- sobol_matrices(N = N, params = params, matrices = matrices)
y <- bratley1992_Fun(mat)
ind <- jansen_ti(d = y, N = N, params = params)

# Plot scatter -------------------------
bratley_scatter <- plot_scatter(data = mat, N = N, Y = y, params = params) +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) + 
  scale_y_continuous(breaks = pretty_breaks(n = 3)) + 
  facet_wrap(~variable, ncol = 6) +
  labs(x = "$x$", y = "$y$") +
  theme(plot.margin = unit(c(0, 0.1, 0, 0), "cm"))

# Plot Ti -----------------------------
bratley_ti <- plot_ti(ind = ind, params = params) + 
  theme(plot.margin = unit(c(0, 0.1, 0, 0), "cm")) + 
  labs(x = "", y = "$T_i$") +
  scale_x_discrete(guide = guide_axis(n.dodge = 1))

# Merge -------------------------------
bratley_plots <- plot_grid(bratley_scatter, bratley_ti, labels = c("f", ""),
                           ncol = 2, rel_widths = c(0.73, 0.27))
bratley_plots


## ----merged_plots_all, dependson=c("bratley_fun", "uncertainty_plots2"), fig.height=5.6, fig.width=6.4, warning=FALSE----

# PLOT BRATLEY ET AL. DISCRETIZED ---------------------------------------------

selected.plots <- 10
bottom.plots <- plot_grid(y.scatters[[selected.plots]] + 
                            facet_wrap(~variable, ncol = 6 ,scales = "free_x") + 
                            scale_x_continuous(breaks = pretty_breaks(n = 3)) + 
                            scale_y_continuous(breaks = pretty_breaks(n = 3))  +
                            labs(x = "$x$", y = "$y$"), 
                          y.indices[[selected.plots]], 
                          ncol = 2, rel_widths = c(0.73, 0.27))

all.bottom.plots <- plot_grid(bratley_plots, bottom.plots, ncol = 1, labels = c("", "g"), 
                              rel_heights = c(0.47, 0.53))
all.bottom.plots


## ----plot_uncertainty, dependson=c("run_model_fun", "arrange_output"), warning = FALSE, fig.height=3, fig.width=3----

# PLOT UNCERTAINTY -----------------------------------------------------------


# SCATTERPLOT ----------------------------------

supp.labs <- c("Symmetric", "Star", "$L_2$", "Centered", "Wrap-around", "Modified")
names(supp.labs) <- discrepancy_methods

# SCATTERPLOT -----------------------------------

a <- melt(final.output, measure.vars = discrepancy_methods) %>%
  .[, variable:= factor(variable, levels =   c("symmetric", "wraparound", 
                                               "centered", "L2",
                                               "star","modified"))] %>%
  ggplot(., aes(cost.discrepancy, k, color = value)) + 
  geom_point(size = 0.4) + 
  scale_colour_gradientn(colours = c("black", "purple", "red", "orange", 
                                     "yellow", "lightgreen"), 
                         name = expression(italic(r)), 
                         breaks = pretty_breaks(n = 3)) +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  labs(x = "Nº of model runs", y = "$k$") + 
  facet_wrap(~variable, ncol = 6, labeller = labeller(variable = supp.labs)) +
  theme_AP() + 
  theme(legend.position = "none", 
        strip.background = element_rect(fill = "white"))

# BOXPLOTS --------------------------------------

discrete.dt <- melt(final.output.discrete, measure.vars = discrepancy_methods) %>%
  .[, Function:= "Bratley et al. 1992"] %>%
  .[, .(variable, value, Function)]

final.dt <- melt(final.output, measure.vars = discrepancy_methods) %>%
  .[, Function:= "Meta-model"] %>%
  .[, .(variable, value, Function)]

b <- rbind(discrete.dt, final.dt) %>%
  ggplot(., aes(reorder(variable, -value), value, fill = Function)) +
  geom_boxplot() + 
  scale_x_discrete(labels = supp.labs) +
  scale_fill_manual(values = wes_palette(2, name = "Chevalier1")) +
  labs(x = "", y = "$r$") + 
  theme_AP() +
  theme(legend.position = "none",
        plot.margin = unit(c(0, 0.1, 0, 0), "cm"))



## ----plot_scatter, dependson="plot_uncertainty", fig.height=4, fig.width=5, warning = FALSE----

# PLOT SCATTERPLOT ------------------------------------------------------------

a



## ----plot_boxplots, dependson="plot_uncertainty", fig.height=2, fig.width=5----

# PLOT BOXPLOTS ---------------------------------------------------------------

b


## ----merged_plots_all2, dependson=c("merge_plots", "bratley_fun", "uncertainty_plots2", "merged_plots_all", "plot_uncertainty"), fig.height=6.6, fig.width=6.4, warning=FALSE, dev="pdf"----

# MERGE ALL PLOTS -------------------------------------------------------------

legend <- get_legend(a + theme(legend.position = "top", 
                               legend.margin=margin(0, 0, 0, 0),
                               legend.box.margin=margin(-7,-7,-7,-7)))

legend.boxplot <- get_legend(b + theme(legend.position = "top", 
                               legend.margin=margin(0, 0, 0, 0),
                               legend.box.margin=margin(-7,-7,-7,-7)))

boxplots.plot <- plot_grid(legend.boxplot, b, ncol = 1, rel_heights = c(0.15, 0.85), 
                           labels = c("e", ""))
scatters.plot <- plot_grid(legend, a, ncol = 1, rel_heights = c(0.18, 0.82), 
                           labels = c("d", ""))
da <- plot_grid(scatters.plot, boxplots.plot, ncol = 1, rel_heights = c(0.55, 0.45))

full.plot <- plot_grid(scat_plot, da, ncol = 2, rel_widths = c(0.25, 0.75))
full.plot

plot_grid(full.plot, all.bottom.plots, ncol = 1, rel_heights = c(0.6, 0.45), 
          labels = c("", ""))



## ----statistics, dependson=c("arrange_output", "arrange_discrete")-------

# STATISTICS ------------------------------------------------------------------

rbind(discrete.dt, final.dt) %>%
  .[, .(mean = mean(value), median = median(value)), .(Function, variable)] %>%
  .[order(-mean, variable)]

















final.output <- fread("final.output.csv")

final.output[, ratio:= cost.jansen / k]
final.output.large <- melt(final.output, measure.vars = discrepancy_methods)

# FOR RANKS-----------
vv <- seq(5, 100, 5)
dt <- lapply(vv, function(x) final.output.large[k <= x, median(value), variable])
names(dt) <- vv

rbindlist(dt, idcol = "k") %>%
  .[, k:= as.numeric(k)] %>%
  ggplot(., aes(k, V1, color = variable)) +
  scale_color_discrete(name = "Estimator") +
  labs(x = expression(italic(k)), 
       y = expression(median(italic(r)))) +
  geom_line() + 
  theme_AP() +
  theme(strip.background = element_rect(fill = "white"))

# Median Nt/k
dt.tmp <- final.output[, .(min = min(ratio), max = max(ratio))]

v <-  seq(0, ceiling(dt.tmp$max), 10)
a <- c(v[1], rep(v[-c(1, length(v))], each = 2), v[length(v)])
indices <- matrix(a, ncol = 2 ,byrow = TRUE)

out.ranks <- list()
for(i in 1:nrow(indices)) {
  out.ranks[[i]] <- final.output[ratio > indices[i, 1] & ratio < indices[i, 2]]
}

names(out.ranks) <- rowMeans(indices)
out.ranks.large <- lapply(out.ranks, function(x) 
  melt(x, measure.vars = discrepancy_methods))

median.mae <- lapply(out.ranks.large, function(x) 
  x[, median(value, na.rm = TRUE), variable]) %>%
  rbindlist(., idcol = "N") %>%
  .[, N:= as.numeric(N)] 


ggplot(median.mae, aes(N, V1, color = variable)) +
  geom_line() + 
  labs(x = "$\frac{C}{d}$", y = "median($r$)") +
  theme_AP() + 
  scale_color_discrete(labels = c("Symmetric", "Star", "$L2$", "Centered", 
                                  "Wrap-around", "Modified")) +
  theme(legend.position = "top")






rbindlist(out.ranks, idcol = "N") %>%
  melt(., measure.vars = discrepancy_methods) %>%
  ggplot(., aes(N, val))