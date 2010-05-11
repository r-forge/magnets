\name{coupling.matrices}
\Rdversion{1.1}
\alias{coupling.matrices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compute the Coupling Matrices between the Different Micro-Spins
}
\description{
  In order to compute the interaction between the different micro-spins,
  the coupling matrices \eqn{\mathcal{C}_{ij}}{C_{ij}} are needed. The interaction energy is

  \deqn{-\frac{1}{2}\sum_{i,j}\bold{m}^{(i)}\cdot\mathcal{C}_{ij}\cdot\bold{m}^{(j)}}{-0.5 * sum(m_i C_{ij} m_j)}

  and the corresponding term in the effective magnetic field
  \eqn{\bold{h}_{eff}^{(i)}}{h_{eff}^(i)} is

  \deqn{\mathcal{C}_{ij}\cdot\bold{m}^{(j)}}{C_{ij} m_j}

This function is used in the function \code{\link{multispin.evol}}, but
can also be used sparately.
}
\usage{
coupling.matrices(r, V)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
    the positions of the micro-spins A \eqn{3\times n}{(3 x n)} array where each
    colume is the spin direction for one spin. \eqn{n} columes for \eqn{n}
    spins. 
}
  \item{V}{
    a \eqn{n} element vector for the volume of each micro-magnet.
}
}
\details{
  The coupling matrices \eqn{\mathcal{C}_{ij}}{C_{ij}} are computed according Eq. (5) in Ref. [1] in \code{\link{multispin.evol}}.
}
\value{
  This function returns the total coupling matrix between the spins with
  dimension \eqn{3m\times 3m}{(3m x 3m)}. It can be regarded as
  \eqn{m\times m}{(m x m)} blocks where each \eqn{3\times 3}{(3 x 3)}
  block is the coupling matrix between any two spins \eqn{i} and
  \eqn{j}. It contains all the information about positions and volumes
  for all spins and external spins. \eqn{m} is the total number of spins
  and external spins.
}
\references{
  check out the references in \code{\link{multispin.evol}}
}
\author{
Hai Qian, Electrical Engineering Department, UCLA
}

\seealso{
  \code{\link{multispin.evol}}  to simulate a micro-spin system and
  generate an evolution as the result and check out the examples there

  \code{\link{free.energy}} to compute the free energy of a micro-spin system
}


