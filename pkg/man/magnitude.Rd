\name{magnitude}
\Rdversion{1.1}
\alias{magnitude}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compute All the Micro-spins' Magnetization Magnitudes
}
\description{
  This function is used to check that the micro-spins keep their
  magnitudes to be exactly 1 during the simulation.
}
\usage{
magnitude(g)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
    an evolution object created by the simulation \code{\link{multispin.evol}}
  }
}
\details{
  Compute \eqn{\vert \bold{m}^{(i)}(t) \vert}{sqrt(m_ix^2 + m_iy^2 +
    m_iz^2)} for all the recorded times.
}
\value{
  The returned result is a "magnitude" object, which is a \eqn{n\times
    l}{n x l} array with each row for one micro-spin in \eqn{l}
  seqential times. \eqn{l} is the length of \code{g$time} defined
  when running the simulation function \code{\link{multispin.evol}}. The
  "magnitude" object also has 3 attributes:
  \item{class}{"magnitude"}
  \item{n}{number of micro-spins, external fixed spins are not included}
  \item{time}{\code{g$time}}
}
\references{
  check out the references in \code{\link{multispin.evol}}
}
\author{
  Hai Qian, Electrical Engineering Department, UCLA
}

\seealso{
  \code{\link{multispin.evol}} to simulate a micro-spin system and
  generate an evolution as the result and check out the examples there
  
  \code{\link{plot.evolution}} to plot various properties of the
  multi-spin system
}


