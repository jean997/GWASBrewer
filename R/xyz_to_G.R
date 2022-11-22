#'@title Generate G from XYZ Specification
#'@param taux_xz,tau_yz Effect size between Z and X or Y as signed percent variance explained (see details)
#'@param dir_xz,dir_yz Effect direction between Z and X or Y (see details)
#'@param beta Signed variance of Y explained by X, see details
#'@details
#'This function generates a matrix G of direct effects corresponding to a model in with variables Y, X, and Z_1, ..., Z_K.
#'There is a causal effect of X on Y given by beta, which specifies the proportion of the variance of Y explained by X.
#'The K variables Z_1, ..., Z_K can have
#'effects to or from X and Y but do not have direct effects on each other. Vectors \code{dir_xz} and \code{dir_yz} specify the direction
#'of these effects with +1 corresponding to "to" effects and -1 corresponding to "from" effects. For example, \code{dir_xz = c(1, -1)}
#'and \code{dir_yz = c(1, 1)} would indicate two variables Z_1 and Z_2 with Z_1 being a common cause of both X and Y and Z_2
#'being a cause of $Y$ and caused by X (Z_2 is a mediator between X and Y). The function will give an error if there
#'is any index with \code{dir_xz} equal to -1 and \code{dir_yz} equal 1. This would indicate a variable that is a mediator between Y and
#'X, however, because there is an effect from X to Y assumed, the resulting graph would be cyclic and not allowed.
#'
#'The inputs `tau_xz` and `tau_yz` specify the effect sizes of each of the Z_k variables on or
#'from X and $Y$. These are given in signed percent variance explained. If we again used \code{dir_xz = c(1, -1)} and \code{dir_yz = c(1, 1)}
#' with \code{tau_xz = c(0.2, -0.3)} and \code{tau_yz = c(0.1, 0.25)},
#' this means that the confounder, Z_1 explains 20\% of the variance of X and 10\% of the
#'variance of Y and both effects are positive. X explains 30\% of the variance of the mediate Z_2 with a negative effect direction and
#'Z_2 explains 25\% of the variance of Y.
#'@return A matrix of direct effects corresponding to variables in the order (Y, X, Z_1, ..., Z_K)
#'@examples
#'xyz_to_G(tau_xz = c(0.2, -0.3), tau_yz = c(0.1, 0.25),
#'         dir_xz = c(1, -1), dir_yz = c(1,1), beta = 0)
#' # code below will give an error due to specification of a cyclic graph.
#'xyz_to_G(tau_xz = c(0.2, -0.3), tau_yz = c(0.1, 0.25),
#'         dir_xz = c(1, -1), dir_yz = c(-1,1), beta = 0.1)
#' # with beta = 0, there is no cycle so no error
#'xyz_to_G(tau_xz = c(0.2, -0.3), tau_yz = c(0.1, 0.25),
#'         dir_xz = c(1, -1), dir_yz = c(-1,1), beta = 0)
#'@export
xyz_to_G <- function(tau_xz, tau_yz, dir_xz, dir_yz, beta){
  nz <- length(tau_xz)
  n <- nz+2
  tau_xz <- check_scalar_or_numeric(tau_xz, "tau_xz", nz)
  tau_yz <- check_scalar_or_numeric(tau_yz, "tau_yz", nz)
  dir_xz <- check_scalar_or_numeric(dir_xz, "dir_nz", nz)
  dir_yz <- check_scalar_or_numeric(dir_yz, "dir_yz", nz)
  if(!all(dir_xz %in% c(-1, 1))){
    stop("dir_xz must have elements -1 or 1 only.")
  }
  if(!all(dir_yz %in% c(-1, 1))){
    stop("dir_yz must have elements -1 or 1 only.")
  }
  tau_xz <- sqrt(abs(tau_xz))*sign(tau_xz)
  tau_yz <- sqrt(abs(tau_yz))*sign(tau_yz)
  beta <- sqrt(abs(beta))*sign(beta)

  # Direct Effects
  G <- matrix(0, nrow = n, ncol = n)
  G[2,1] <- beta
  G[which(dir_xz == 1) + 2, 2] <- tau_xz[dir_xz ==1]
  G[2, which(dir_xz == -1) + 2] <- tau_xz[dir_xz == -1]

  G[which(dir_yz == 1) + 2, 1] <- tau_yz[dir_yz ==1]
  G[1, which(dir_yz == -1) + 2] <- tau_yz[dir_yz == -1]
  Gc <- check_G(G, rep(1, n), n)
  return(G)
}
