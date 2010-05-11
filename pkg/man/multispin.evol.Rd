\name{multispin.evol}
\Rdversion{1.1}
\alias{multispin.evol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Numerical Simulation of 3D Micro-magnets in 3D Dimensions under A Finite Temperature
}
\description{
  Numerical simulates the time evolution of the rescaled
  Landau-Lifshitz-Gilbert (LLG) equation for any number of 3D
  micro-magnets in 3D dimensional space under a finite temperature.
  External magnetic field and
  external fixed micro-magnets are allowed. The magnitude of the spins
  is always 1 because of the rescaling of the LLG equation.
}
\usage{
multispin.evol(init = NULL, positions = NULL, volumes = NULL, fixed.ext
= NULL, pos.ext = NULL, vol.ext = NULL, coupling.matrix = NULL, T = 10000, dt = 0.01, extern.field = 0, biaxial.anisotropy = 0,
demag.tensor = 0, given.noise = TRUE, noise.type = "gaussian", sd = 1,
mean = 0, dw = NULL, record = TRUE, record.step = 1, record.list = NULL,
alpha = 0.1, nu = 0.02, error = 10^-10, iter = 1000, implicit = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{init}{
    the initial configuration of the micro-spins. A \eqn{3\times n}{(3 x n)} array where each
    colume is the spin direction for one spin. \eqn{n} columes for \eqn{n}
    spins. Must be given for the simulation.
  }
  
  \item{positions}{
    the positions of the spins. A \eqn{n\times 3}{(n x 3)} array where each row is the \eqn{xyz}
    position for one spin. The order for the spins is the same with
    that in init. Can be NULL if \code{coupling.matrix} is given.
  }
  
  \item{volumes}{
    a \eqn{n} element vector for the volume of each micro-magnet. Can be NULL
    if \code{coupling.matrix} is given.
  }
  
  \item{fixed.ext}{
    the configuration of the external micro-magnets. It takes the same
    data structure as init, but will be absolutely fixed in the
    simulation. Can be NULL if there is no external fixed magnets in the system.
  }
  
  \item{pos.ext}{
    the positions of the fixed external micro-magnets similar to
    positions. Can be NULL if \code{coupling.matrix} is given.
  }
  
  \item{vol.ext}{
    the volumes of the fixed external micro-magnets similar to
    volumes. Can be NULL if \code{coupling.matrix} is given.
  }
  
  \item{coupling.matrix}{
    the total coupling matrix between the spins with dimension
    \eqn{3m\times 3m}{(3m x 3m)}. It can be regarded as \eqn{m\times
    m}{(m x m)} blocks where each \eqn{3\times 3}{(3 x 3)} block
    is the coupling matrix between any two spins \eqn{i} and \eqn{j}. It contains all the
    information about positions and volumes for all spins and external
    spins, and this is why \code{positions} and \code{volumes} can be NULL when this
    quantity is given. \eqn{m} is the total number of spins and external
    spins. If \code{coupling.matrix}, \code{positions} and \code{volumes} are all given, the
    latter two will be ignored. If \code{coupling.matrix} is NULL, then \code{positions} and
    \code{volumes} must be given. \code{coupling.matrix} is computed with the function \code{\link{coupling.matrices}}
  }
  
  \item{T}{
    an integer, the total number of time steps to evolve the system. 
  }
  
  \item{dt}{
    time step size.
  }
  
  \item{extern.field}{
    the external magnetic field. Can be a 3-vector representing \eqn{xyz}
    components of the external magnetic field, if it is uniform and
    time-independent. Can be a \eqn{3n}{3n} vector with every 3 elements for each
    spin, if it is non-uniform and time-independent. Can be a \eqn{3T}{3T} vector
    if it is uniform and time-dependent. Can be a \eqn{3nT}{3nT} vector if it is
    non-uniform and time-dependent.
  }
  
  \item{biaxial.anisotropy}{
    a \eqn{n} vector with each element being the biaxial anistropy for each
    spin. A single number can be used here if all spins have the same
    biaxial anisotropy.
  }

  \item{demag.tensor}{
    the effective demagnetization tensor \eqn{N_x}{Nx}, \eqn{N_y}{Ny}, \eqn{N_z}{Nz}. The uniaxial
    anisotropy should be already included in the effective
    demagnetization tensor. Can be a single number if all components for
    all spins are the same; can be a 3-vector if all spins are the same;
    can be a \eqn{3n}{3n} vector where every 3-element corresponds to each spin.
  }
  
  \item{given.noise}{
    whether to use external noise for test purpose, or use faster
    internal noise generation.
  }
  
  \item{noise.type}{
    noise type. For now, only Gaussian noise is implemented. More are
    coming in the future.
  }
  
  \item{sd}{
    noise standard deviation. obsolete argument now. Possibly use for
    future use for different type of noise.
  }
  
  \item{mean}{
    noise mean value. obsolete argument now. Possibly use for
    future use for different type of noise.
  }
  
  \item{dw}{
    a \eqn{3n\times T}{(3n x T)} array. Give noise for test purpose. If
    it is NULL, the noise will be generated internally in the routine.
  }
  
  \item{record}{
    whether to record the spin configurations.
  }
  
  \item{record.step}{
    every \code{record.step} to record the spin configurations. This takes
    effect only when \code{record.list} is NULL. The default value is 1 and the
    configurations at every time step are recorded.
  }

  \item{record.list}{
    an integer vector whose elements are the time step to record the
    spin configurations. If it is NULL, the record times are decided by
    \code{record.step}. 
  }
  
  \item{alpha}{
    the dissipation coefficient. The typical value of \eqn{0.001-0.1}.
  }
  
  \item{nu}{
    The intensity of the thermal noise. Decided by the temperature.
  }
  
  \item{error}{
    in the implicit method, the error requirement for solving the
    implicit equations at each time step.
  }
  
  \item{iter}{
    in the implicit method, the maximum iteration number for solving teh
    implicit equations at each time step.
  }
  
  \item{implicit}{
    whether to use implicit or explicit method. The implicit method is more
    accurate and the error can be controlled, while the explicit method is
    much much faster. See the detailed explanation in the next section.
  }
}
\details{
  \strong{(It is much better to view this section in the PostScript version
  of this documentation. To do that, just run "help (multispin.evol,
    help_type="pdf")" in R command line after loading this library and you will find "multispin.evol.pdf" in the current working directory.)}

    \emph{1. The Simulation Equations}
    
  With all the necessary inputs in the above section, this function uses
    the LLG equation to drive the system to evolve. The numerical method
    employed here is basically the mid-point method, which is able to
    keep the magnitude of the spins to be exactly 1 in the whole
    process.

    The original rescaled LLG equation with thermal noise is

    \deqn{\frac{d\bold{m}^{(i)}}{dt}=-\bold{m}^{(i)}\times
    \bold{h}_{eff}^{(i)}\left(\left\{\bold{m}^{(j)}\right\}\right) +
    \alpha\bold{m}^{(i)}\times\frac{d\bold{m}^{(i)}}{dt}-\bold{m}^{(i)}\times\nu
    d\bold{w}^{(i)}\qquad\qquad(1)}{dm_i/dt = -m_i x h_{eff}({m_j}) +
    alpha m_i x dm_i/dt - m_i x nu dw_i     (1)}  

    where \eqn{\bold{h}_{eff}^{(i)}}{h_{eff}} depends on all the
    micro-spins in the system and is given by

    \deqn{\bold{h}_{eff}^{(i)}=\bold{h}_{ext} - \mathcal{N}_{eff}^{(i)}\cdot
      \bold{m}^{(i)} + \sum_i \mathcal{C}_{ij}\cdot\bold{m}^{(j)} -
      2k_2\sum_{\alpha=x,y,z}m^{(i)}_{\alpha}\left(1-(m^{(i)}_{\alpha})^2\right)\hat{\bold{e}}_{\alpha}\qquad\qquad(2)}{h_{eff}=h_{ext}-N_{eff}
      m_i + sum_i C_ij m_j - 2 k_2 sum(m_a(1-m_a^2))        (2)}

    and \eqn{d\bold{w}^{(i)}}{dw_i} is the noise with standard deviation
    equal to 1. The above equation Eq. (1) is equivalent to

    \deqn{\frac{d\bold{m}^{(i)}}{dt}=-\frac{1}{1+\alpha^2}\bold{m}^{(i)}\times\left(
    \bold{h}_{eff}^{(i)}+\nu d\bold{w}^{(i)}\right)
    -\frac{\alpha}{1+\alpha^2}\bold{m}^{(i)}\times\left[\bold{m}^{(i)}\times\left(\bold{h}_{eff}^{(i)}+\nu d\bold{w}^{(i)}
    \right)\right]\qquad(3)}{dm_i/dt = -(1+alpha^2)^(-1) m_i x
    (h_{eff}+nu dw_i) -alhpa/(1+alpha^2) m_i x [m_i x (h_{eff} + nu
    dw_i)]        (3)}

    In the so-called mid-point method, the LLG equation Eq. (1) is

    \deqn{\frac{\bold{m}^{(i)}(t+dt)-\bold{m}^{(i)}(t)}{dt} =
      -\frac{\bold{m}^{(i)}(t+dt)+\bold{m}^{(i)}(t)}{2}\times
      \bold{h}_{eff}^{(i)}\left(\left\{\frac{\bold{m}^{(j)}(t+dt)+\bold{m}^{(j)}(t)}{2}\right\}\right)+   
      }{(m_i(t+dt)-m_i(t))/dt = (m_i(t+dt)+m_i(t))/2 x
	h_{eff}({(m_j(t+dt)+m_j(t))/2}) +}
      \deqn{\quad\quad\quad\quad\quad +
    \alpha\frac{\bold{m}^{(i)}(t+dt)+\bold{m}^{(i)}(t)}{2}\times\frac{\bold{m}^{(i)}(t+dt)-\bold{m}^{(i)}(t)}{dt}
    - }{
    + alpha (m_i(t+dt)+m_i(t))/2 x (m_i(t+dt)-m_i(t))/dt - }
  \deqn{ -
    \frac{\bold{m}^{(i)}(t+dt)+\bold{m}^{(i)}(t)}{2}\times\nu d\bold{w}^{(i)}\qquad\qquad(4)}{-
    (m_i(t+dt)+m_i(t))/2 x nu dw_i        (4)}

     Then we have to solve a set of \eqn{3n}{3n} entangled nonlinear equation in one
    shot. This is the implicit method. The computation is quite demanding.

    The other version of the LLG equation Eq. (3) can be re-write to remove the time-derivative
    term as

    \deqn{\frac{\bold{m}^{(i)}(t+dt)-\bold{m}^{(i)}(t)}{dt} =
    -\frac{1}{1+\alpha^2}\frac{\bold{m}^{(i)}(t+dt)+\bold{m}^{(i)}(t)}{2}\times\left(\bold{h}_{eff}^{(i)}(\cdot)+\nu
    d\bold{w}^{(i)}\right)
    - }{(m_i(t+dt)-m_i(t))/dt = (1+alpha^2)^(-1) (m_i(t+dt)+m_i(t))/2 x
    (h_{eff}(.) + nu dw_i) - }
  \deqn{\quad\quad\quad\quad\quad -
    \frac{\alpha}{1+\alpha^2}\frac{\bold{m}^{(i)}(t+dt)+\bold{m}^{(i)}(t)}{2}\times
    \left[\bold{m}^{(i)}(t)\times\left(\bold{h}_{eff}^{(i)}(\cdot)+\nu
    d\bold{w}^{(i)}\right)\right]\qquad(5)}{- alpha/(1+alpha^2)
    (m_i(t+dt)+m_i(t))/2 x [m_i(t) x (h_{eff}(.) + nu)]       (5)}
    
    If we replace the arguments of the effective magnetic field \eqn{\bold{h}_{eff}^{(i)}}{h_{eff}}
    with just \eqn{\{m_j(t)\}}{{m_j(t)}}, then we just have to solve \eqn{n}{n} sets of equations
    with each set having 3 linear equations. This is the explicit method.
    And it reduced the computation amount and saved time. As long as the
    time step is small enough, it should be the same with the implicit
    method. 

    By default, this function uses explicit method.

    With a given noise, one can compare the simulation result of these
    two method under a certain temperature.

    \emph{2. The Rescaled Quantities}
    
    In the above equations, the quantities are all rescaled

    \deqn{\frac{\bold{M}^{(i)}}{M_s} \rightarrow \bold{m}^{(i)}}{M_i /
    Ms -> m_i}
    \deqn{\frac{\bold{H}_{ext}^{(i)}}{M_s} \rightarrow
    \bold{h}_{ext}^{(i)}}{h_{ext} = H_{ext} / M_s}
  \deqn{t\rightarrow \gamma M_s t}{t -> gamma Ms t}
  \deqn{\frac{K_2}{\mu_0 M_s^2}\rightarrow k_2}{K2/(mu_0 Ms^2) -> k2}
    \deqn{\nu = \sqrt{\frac{2\alpha k_B T}{\mu_0 M_s^2 V}}}{nu = sqrt(2
      alpha k_B T/(mu_0 Ms^2 V))}
    \deqn{\frac{F}{\mu_0 M_s^2 V}\rightarrow f}{F/(mu_0 Ms^2 V) -> f}

    where \eqn{V} is the unit volume and \eqn{M_s}{Ms} is the saturation
    magnetization.

    \emph{3. The SI Units of Some Quantities}

    \eqn{[\mu_0]=N/A^2}{mu_0 = N/A^2}, \eqn{[\gamma]=m/(A\cdot
      s)}{[gamma] = m/(A s)}, \eqn{[Ms]=A/m}{[Ms] = A/m},
      \eqn{[K_1]=[K_2]=J/m^3}{[K1] = [K2] =J/m^3}
}
\value{
  The output of this function is an evolution class object, which is a
  list with two components: process and time:

  \item{process}{
    a \eqn{3n\times\mbox{length(record.list)}}{(3n x
    length(record.list))} array describing the \eqn{3n}{3n}
    components of all \eqn{n}{n} spins.
  }
  \item{time}{
    the same as record.list+1, an integer vector. The time starts from 1.
  }

  \item{}{This evolution class object also has several attributes:}

  \item{class}{"evolution"}
  \item{call}{ explicit or implicit method }
  \item{dt}{ the time step size }
  \item{n}{ the number of spins in the system }
  \item{vol}{ volumes of the \eqn{n} spins }
  \item{pos}{ positions of the \eqn{n} spins }
  \item{ext.n}{ external fixed spins }
  \item{ext.spins}{ external spin configuration }
  \item{ext.vol}{ external spin volumes }
  \item{ext.pos}{ external spin positions }
  \item{K2}{ the biaxial anisotropy parameters }
  \item{eff.N}{ the effective demagnetization tensor }
  \item{coupling}{ the coupling matrix including all internal and
  external spins }
}
\references{
  [1] \emph{Performance of Magnetic Quantum Cellular Automata and Limitations due
  to Thermal Noise}, F. M. Spedalieri, A. P. Jacob, D. Nikonov, and V. P.
  Roychowdhury, arXiv: 0906.5172, \url{http://arxiv.org/abs/0906.5172}

  [2] \emph{Nonlinear Magnetization Dynamics in Nanosystems}, G. Bertotti,
  I. Mayergoyz andC. Serpico, 2009 Elsevier Ltd. (Chap. 11 Numerical
  techniques for magnetization dynamics analysis)

  [3] \emph{Nonlinear Magnetization Dynamics in Thin-films and
    Nanoparticles}, Massimiliano d'Aquino, Doctorate thesis in
    Electrical Engineering, University of Napoli, 2004. 
}
\author{
Hai Qian, Electrical Engineering Department, UCLA
}
\note{
  All the simulation is done in orthogonal coordinates \eqn{xyz}, and
    the easy axis is along one of them. Thus the demagnetization tensor
    is simplified to a \eqn{3\times 3}{(3 x 3)} diagonal matrix where
    all non-diagonal elements are 0.
  
  The time evolution starts from 0, i.e. the smallest value in record.list
  can be 0.

  One must keep the same order for the micro-spins when defining
  \code{init}, \code{positions}, \code{volumes}, \code{coupling.matrix}
  ... etc.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{coupling.matrices}} to compute the coupling matrix between
  the spins

  \code{\link{rect.demag}} to compute the demanetization
  tensor of a rectangular prism given its 3-dimension sizes
  
  \code{\link{magnitude}} to check whether the magnitude of the spins is
  exactly 1
  
  \code{\link{free.energy}} to compute the free energy of the system
  
  \code{\link{component.animation}} to display an animation for one
  component of all spins

  \code{\link{plot.evolution}} to plot various properties of an
  evolution object
}
\examples{
#
# Example 1: simulate a string of 25 micro-spins
#
a <- 10; b <- 36; c <- 6; # unit is nm
V0 <- a*b*c * (10^-9)^3 # m^3

kT <- 0.026 * 1.61 * 10^(-19); # J
Ms <- 796 * 10^3; # A/m
gamma <- 2.21 * 10^5; # m/(A s)
mu0 <- 4*pi*10^(-7); # N/A^2
alpha <- 0.05; # dissipation coefficient
K1 <- 3.2 * 10^3; # J/m^3

N.v0 <- c(rect.demag(b,c,a), rect.demag(c,a,b), rect.demag(a,b,c));
N.v <- N.v0;
N.v[2] <- N.v0[2] - 2*K1/(mu0*Ms^2); # effective demagnetization tensor

nu <- sqrt (2*alpha*kT/(mu0*Ms^2*V0));
k2.min <- 0.5 * (N.v[1] - N.v[2]);

n <- 25; # number of micro-magnets

d <- 15;
x <- seq (from=0, by=d, length.out=n);
y <- rep (0, n);
z <- rep(0, n);
r <- cbind (x,y,z); # positions of n spins

V <- rep (a*b*c, n); # volumes of the magnets
m <- array (rep (c(1,0,0),n), dim=c(3,n)); # initial configuration

xe <- c(-d);
ye <- c(0);
ze <- c(0);
re <- cbind (xe, ye, ze); # external magnet positions
me <- array (c(0,1,0), dim=c(3,1)); # external magnet configuration
v.ext <- c(2200);

ts <- gamma * Ms; # time scale
T <- 40000; # 1 ps is equal to 10^(-12)*ts
dt <- 0.01;

h.ext <- rep (0, 3*n); # external field

pos <- rbind (r, re); # external magnets
vol <- c(V, v.ext);
Cp <- coupling.matrices (pos, vol);

g <- multispin.evol (init=m, coupling.matrix=Cp, T=T, dt=dt,
fixed.ext=me, extern.field=h.ext, demag.tensor=N.v, nu=0.01395,
biaxial.anisotropy=0.3*k2.min, alpha=0.05, given.noise=FALSE,
implicit=FALSE);

fg <- free.energy (g); # free energy
plot (fg);

mg <- magnitude (g);
plot (g, spins=c(1,2,5,10)) # plot magnitudes
plot (g, spins=c(1,2,5,10), components=c(2,2,2,3)); # plot time-dependence of several spins

ani <- component.animation (g, by=10, time=T, component=2); # plot animations
replot.animation (ani);

# ------------------------------------------
#
# Example 2: simulate a 5 spin gate
#
a <- 15; b <- 30; c <- 6; # nm
V0 <- a*b*c * (10^-9)^3 # m^3

kT <- 0.026 * 1.61 * 10^(-19); # J
Ms <- 796 * 10^3; # A/m
gamma <- 2.21 * 10^5; # m/(A s)
mu0 <- 4*pi*10^(-7); # N/A^2
alpha <- 0.05;
K1 <- 3.2 * 10^3; # J/m^3

N.v0 <- c(rect.demag(b,c,a), rect.demag(c,a,b), rect.demag(a,b,c));
N.v <- N.v0;
N.v[2] <- N.v0[2] - 2*K1/(mu0*Ms^2);

nu <- sqrt (2*alpha*kT/(mu0*Ms^2*V0));
k2.min <- 0.5 * (N.v[1] - N.v[2]);

n <- 5;

# positions of 5 spins
d <- 25;
x <- c(0, -d, 0, 0, d);
y <- c(d*2^(1/3), 0, -d*2^(1/3), 0, 0);
z <- rep(0, 5);
r <- cbind (x,y,z);

V <- rep (a*b*c, n); # volumes of the magnets
m <- array (rep (c(1,0,0),n), dim=c(3,n)); # initial configuration

# positions of external spins, volumes can be different
# 3 external spins totally
xe <- c(0, -2*d, 0);
ye <- c(2*d*2^(1/3), 0, -2*d*2^(1/3));
ze <- c(0,0,0);
re <- cbind (xe, ye, ze);

me <- cbind (c(0,1,0), c(0,1,0), c(0,-1,0));
v.ext <- rep (3000, 3);

ts <- gamma * Ms;
T <- 35183; # 1 ps is equal to 10^(-12)*ts
dt <- 0.01;

# time-dependent external field
tmp <- seq (0.15, 0, length.out=30000);
x.i <- c(tmp, rep(0,T-30000));
h.ext <- c();
for (i in 1:T) {
  h.ext <- c(h.ext, rep(c(x.i[i],0,0), n));
}

pos <- rbind (r, re);
vol <- c(V, v.ext);
Cp <- coupling.matrices (pos, vol);

g <- multispin.evol (init=m, coupling.matrix=Cp, T=T, dt=dt,
fixed.ext=me, extern.field=h.ext, demag.tensor=N.v, nu=0.01395,
biaxial.anisotropy=0.0, alpha=0.05, given.noise=FALSE, implicit=FALSE);

fg <- free.energy (g);
plot (fg);

mg <- magnitude (g);

ani <- component.animation (g, by=10, time=T, component=2);
replot.animation (ani);

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~interface }

