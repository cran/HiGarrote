% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HiGarrote.R
\name{HiGarrote}
\alias{HiGarrote}
\title{An Automatic Method for the Analysis of Experiments using Hierarchical Garrote}
\usage{
HiGarrote(
  D,
  y,
  quali_id = NULL,
  quanti_id = NULL,
  heredity = "weak",
  U = NULL,
  me_num = NULL,
  quali_contr = NULL
)
}
\arguments{
\item{D}{An \eqn{n \times p} data frame for the unreplicated design matrix, where \eqn{n} is the run size and \eqn{p} is the number of factors.}

\item{y}{A vector for the responses corresponding to \code{D}. For replicated experiments, \code{y} should be an \eqn{n \times r} matrix, where \eqn{r} is the number of replicates.}

\item{quali_id}{A vector indexing qualitative factors.}

\item{quanti_id}{A vector indexing quantitative factors.}

\item{heredity}{Choice of heredity principles: weak or strong. The default is weak.}

\item{U}{Optional. An \eqn{n \times P} model matrix, where \eqn{P} is the number of potential effects. 
The inclusion of potential effects supports only main effects and two-factor interactions. 
Three-factor and higher order interactions are not supported.
The colon symbol ":" must be included in the names of a two-factor interaction for separating its parent main effects.
By default, \code{U} will be automatically constructed.
The potential effects will then include all the main effects of qualitative factors, the first two main effects (linear and quadratic) of all the quantitative factors, and all the two-factor interactions generated by those main effects.
By default, the coding systems of qualitative and quantitative factors are Helmert coding and orthogonal polynomial coding, respectively.}

\item{me_num}{Optional. A \eqn{p \times 1} vector for the main effects number of each factor.
\code{me_num} is required when \code{U} is not \code{NULL} and must be consistent with the main effects number specified in \code{U}.}

\item{quali_contr}{Optional. A list specifying the contrasts of factors. 
\code{quali_contr} is required only when the main effects of a qualitative factor are not generated by the default Helmert coding.}
}
\value{
A vector for the nonnegative garrote estimates of the identified effects.
}
\description{
`HiGarrote()` provides an automatic method for analyzing experimental data. 
This function applies the nonnegative garrote method to select important effects while preserving their hierarchical structures.
It first estimates regression parameters using generalized ridge regression, where the ridge parameters are derived from a Gaussian process prior placed on the input-output relationship. 
Subsequently, the initial estimates will be used in the nonnegative garrote for effects selection.
}
\examples{
# Cast fatigue experiment
data(cast_fatigue)
X <- cast_fatigue[,1:7]
y <- cast_fatigue[,8]
HiGarrote::HiGarrote(X, y)

# Blood glucose experiment
data(blood_glucose)
X <- blood_glucose[,1:8]
y <- blood_glucose[,9]
HiGarrote::HiGarrote(X, y, quanti_id = 2:8) 

\donttest{
# Router bit experiment
data(router_bit)
X <- router_bit[, 1:9]
y <- router_bit[,10]
for(i in c(4,5)){
my.contrasts <- matrix(c(-1,-1,1,1,1,-1,-1,1,-1,1,-1,1), ncol = 3)
X[,i] <- as.factor(X[,i])
contrasts(X[,i]) <- my.contrasts
colnames(contrasts(X[,i])) <- paste0(".",1:(4-1))
}
U <- model.matrix(~.^2, X)
U <- U[, -1]  # remove the unnecessary intercept terms from the model matrix
me_num = c(rep(1,3), rep(3,2), rep(1, 4))
quali_contr <- list(NULL, NULL, NULL,
                    matrix(c(-1,-1,1,1,1,-1,-1,1,-1,1,-1,1), ncol = 3),
                    matrix(c(-1,-1,1,1,1,-1,-1,1,-1,1,-1,1), ncol = 3),
                    NULL, NULL, NULL, NULL)
HiGarrote::HiGarrote(X, y, quali_id = c(4,5), U = U, 
me_num = me_num, quali_contr = quali_contr)

# Experiments with replicates
# Generate simulated data
data(cast_fatigue)
X <- cast_fatigue[,1:7]
U <- data.frame(model.matrix(~.^2, X)[,-1])
error <- matrix(rnorm(24), ncol = 2) # two replicates for each run
y <- 20*U$A + 10*U$A.B + 5*U$A.C + error
HiGarrote::HiGarrote(X, y)
}

}
\references{
Yu, W. Y. and Joseph, V. R. (2024) "Automated Analysis of Experiments using Hierarchical Garrote," arXiv preprint arXiv:2411.01383.
}
