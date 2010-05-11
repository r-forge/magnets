\name{rect.demag}
\Rdversion{1.1}
\alias{rect.demag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compute the Demagnetization Tensor for A Rectanular Prism
}
\description{
  This function computes the demagnetization tensor for a rectangular
  prism in orthogonal coordinates \eqn{xyz}. The easy axis of the prism
  is along one the coordinate axis. Under such settings, the
  demagnetization tensor is a \eqn{3\times 3}{(3 x 3)} diagonal
  matrix. And this function computes the 3 diagonal non-zero elements
  \eqn{N_x}{Nx}, \eqn{N_y}{Ny} and \eqn{N_z}{Nz}.
}
\usage{
rect.demag(a, b, c)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
    the size along the direction which you do not want to compute the
    demagnetization tensor.
}
  \item{b}{
    this is symmetric to \code{a}
}
  \item{c}{
    the size along the direction which you want to compute the
    demagnetization tensor. For example, to compute \eqn{N_z}{Nz} one
    just puts \eqn{L_z}{Lz} of the prism as \code{c}.
}
}

\value{
  This function returns one of the three non-zero
  elements. \eqn{N_x}{Nx} if \code{c}=\eqn{L_x}{Lx}; \eqn{N_y}{Ny} if
  \code{c}=\eqn{L_y}{Ly}; \eqn{N_z}{Nz} if \code{c}=\eqn{L_z}{Lz}.
}
\references{
  check out the references in \code{\link{multispin.evol}}
}
\author{
  Hai Qian, Electrical Engineering Department, UCLA
}

\seealso{
  \code{\link{multispin.evol}} uses demagnetizatoin tensor which can be
  assigned by hand or calculated with this function given the shape and
  sizes of the micro-magnets. Check out the examples there.
}


