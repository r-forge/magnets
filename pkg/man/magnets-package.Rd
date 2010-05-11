\name{magnets-package}
\Rdversion{1.1}
\alias{magnets-package}
\alias{magnets}
\docType{package}
\title{
  Numerically Simulate n Micro-spins under the Influence of Temperature &
  Other External Conditions
}
\description{
  In recent years, micro-magnets systems have arisen much theoretical and experimental
interest as a possible candidate of the next-generation circuits to
replace the current electronic
integrated circuits. The micro-magnets on the sizes of micro-meters can be
arranged to form logic gates and even wires to transfer
information. Although the speed of micro-magnet components is 
slower than the traditional electronic circuit, the energy that is
needed to change the state of a micro-magnet is much smaller. This
makes it suitable for many special cases where low energy consumption
is the priority. An ideal integrated circuit would be a mix of
the micro-magnet circuit and the traditional electronic
circuit. 

Therefore it is important to be able to do numerical simulations about
any configurations of micro-magnets. However there is no package in R
that can do this job yet. My package will fulfill this gap. 

As the first step, my package now is able to simulate the
time-dependent dynamics of n magnet spins with any possible
configurations. My package has the ability to display the dynamics as
an animation with the help of rggobi package, which proves to be very
helpful in my research. Some system related
properties such as the free energy and the average magnitude of the spins can be
measured.

So far this package can be used as a first step simulation to study
the dynamics of micro-magnet wire (a linear array of micro-magnets) and logic gates. The users are
required to input the configuration of the system, i.e. how the spins
are positioned in the 3D space. Different configurations can have
different dynamics, and thus may be used as different functional
components in an integrated circuit.

In the future versions, More functions will be added. And the package will be able to simulate the detailed dynamics
of one single magnet spin.

This package should be useful to physicists and electrical engineers,
who are interested in design the next generation integrated circuits.
  
  Numerically simulate \eqn{n} micro-spins under certain temperature,
  external magnetic field, external fixed magnets. The
  Landau-Lifshitz-Gilbert Equation is used to simulate the system.
}
\details{
\tabular{ll}{
Package: \tab magnets\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2009-10-15\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}

Use \code{\link{multispin.evol}} to simulate, use
\code{\link{free.energy}} to calculate the free energy of the system,
use \code{\link{component.animation}} and \code{\link{replot.animation}}
to view the dynamics directly as an animation.

}
\author{
  Hai Qian, Electrical Engineering Department, UCLA

Maintainer: Hai Qian <haiqian@ee.ucla.edu>
}
\references{
  check out the references in \code{\link{multispin.evol}}
}
\keyword{ package }

