\name{single.spin.energy}
\Rdversion{1.1}
\alias{single.spin.energy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compute the Free Energy of A Single Isolated Micro-spin
}
\description{
  Calculate the free energy of a single isolated micro-spin
}
\usage{
single.spin.energy(m, N.v, h.ext, K2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
    a 3-vector for spin configuration
  }
  \item{N.v}{
    a 3-vector for effective demagnetization tensor, by properly choosing
    coordinates system to make the tensor diagonal and have only 3
    non-zero elements. The uniaxial anisotropopy should be already included.
  }
  \item{h.ext}{
    3-vector for the external magnetic field
  }
  \item{K2}{
    the biaxial anisotropy parameter
  }
}
\details{
  Please refer to \code{\link{free.energy}} by ignoring the interaction energy
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

  \code{\link{spin.energy.config}} to plot a single isolated spin free
  energy v.s. angle away from \eqn{x}-axis
}


