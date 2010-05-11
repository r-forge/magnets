\name{spin.energy.config}
\Rdversion{1.1}
\alias{spin.energy.config}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot A Single Isolated Spin Free Energy v.s. Angle to x-axis
}
\description{
  This function plot the free energy of a single isolated spin v.s. the
  angle to \eqn{x}-axis. The spin is restricted in \eqn{xy}-plane.
}
\usage{
spin.energy.config(h.ext = c(0, 0, 0), K2 = 0, N.v = c(0, 0, 0), range = seq(-1.3 * pi, 1.3 * pi, 0.02))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h.ext}{
    a 3-vector for external magnetic field
}
  \item{K2}{
    the biaxial anisotropy parameter
}
  \item{N.v}{
    a 3-vector for effective demagnetization tensor, by properly choosing
    coordinates system to make the tensor diagonal and have only 3
    non-zero elements. The uniaxial anisotropopy should be already included.
}
  \item{range}{
    the range of angle 
}
}

\references{
  check out the references in \code{\link{multispin.evol}} 
}
\author{
  Hai Qian, Electrical Engineering Department, UCLA
}

\seealso{
  \code{\link{multispin.evol}} for examples
  
  \code{\link{free.energy}} to compute multiple spin system free energy

  \code{\link{single.spin.energy}} is employed by this function to
  compute the free energy for a single isolated spin
}
