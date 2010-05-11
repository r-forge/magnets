# calculate the properties of the multi-spin system

# input is a evolution object
# ext whether to include external magnets
# when ext is TRUE, include the interaction energy between fixed and physical magnets
# ext.field, whether to include the interaction between external magnets and the external field
free.energy <- function (g, h.ext=NULL, ext=FALSE, ext.field=FALSE)
  {
    l <- length (g$time);
    n <- attr (g, "n");
    n3 <- 3 * n;
    ne <- attr (g, "ext.n");
    fg <- array (0, dim=c(l,2));

    N.v <- attr (g, "eff.N");
    #h.ext <- attr (g, "h.ext");
    Cp <- attr (g, "coupling");
    K2 <- rep (attr (g, "K2"), each=3);
    me <- as.vector (attr (g, "ext.spins"));
    dt <- attr (g, "dt");

    if (is.null(h.ext)) {
      flag <- TRUE;
      h.ext <- rep (0, n3);
    } else if (length(h.ext) == n3) {
      flag <- TRUE;
    } else if (length(h.ext) == n3 * l) {
      flag <- FALSE;
    }
    
    fg[,1] <- g$time*dt;
    for (i in 1:l) {
      if (flag) {
        tmp <- h.ext;
      } else {
        tmp <- h.ext[((g$time[i]-1)*n3+1):(g$time[i]*n3)];
      }
      cmp <- Cp %*% c(g$process[,i],me);
      h <- tmp - 0.5 * N.v * g$process[,i] + 0.5 * cmp[1:n3] - K2 * g$process[,i] * (1-0.5*g$process[,i]^2);
      fg[i,2] <- - sum (g$process[,i]*h);
      if (ext) {
        fg[i,2] <- fg[i,2] - sum(0.5 * me * cmp[(n3+1):(n3+3*ne)]);
      }
      if (ext.field) { # suppose the field is same everywhere
        h <- rep (tmp[1:3], ne);
        fg[i,2] <- fg[i,2] - sum(h * me); 
      }
    }

    class (fg) <- "free.energy";
    return (fg);
  }

# ------------------------------------------------------------------------
single.spin.energy <- function (m, N.v, h.ext, K2)
  {
    h <- h.ext - 0.5 * N.v * m - K2 * m * (1-0.5*m^2);
    return (-sum(m*h));
  }

# ----------------------------------------------------------------------
# check magnitude
magnitude <- function (g)
  {
    l <- length (g$time);
    mg <- array (0, dim=c(n,l));

    class (mg) <- "magnitude";
    attr (mg, "time") <- g$time;
    n <- attr (g, "n");
    
    for (i in 1:n) {
      i1 <- (i-1)*3 + 1;
      i2 <- i1 + 1;
      i3 <- i2 + 1;
      mg[i,] <- sqrt (g$process[i1,]^2 + g$process[i2,]^2 + g$process[i3,]^2);
    }

    return (mg);
  }

# ----------------------------------------------------------------------
summary.evolution <- function (object, ...)
  {
    att <- attributes (object);
    rst <- list (n=att$n, time=max(object$time), dt=att$dt, call=att$call);

    class (rst) <- "summary.evolution";
    rst
  }

# ----------------------------------------------------------------------
print.summary.evolution <- function (x, ...)
  {
    cat ("Call: ");
    print (x$call);
    cat ("\n");

    cat ("number of spins: ");
    print (x$n);
    cat ("\n");

    cat ("evolution time: ");
    print (x$time);
    cat ("\n");

    cat ("time step: ");
    print (x$dt);
    cat ("\n");
  }

print.evolution <- function (x, ...)
  {
    summary.evolution (x)
  }

# ---------------------------------------------------------------------
plot.free.energy <- function (x, ...)
  {
    plot.default (x, type='l', xlab="Time", ylab="Free Energy", main="The Free Energy of the Multi-Spin System", ...);
  }

# --------------------------------------------------------------------
plot.evolution <- function (x=NA, spins=NULL, components=NULL, ...)
  {
    if (is.null(spins) && is.null(components)) {
      fg <- free.energy (x);
      plot (fg);
      
    } else if (is.null(components)) {
      mg <- magnitude (x);
      plot (x$time, mg[spins[1],], ylim=c(0,2), type='l', xlab="Time", ylab="Magnitude", main="Magnitudes of the Spins");
      if (length (spins) > 1) {
        for (i in 2:length(spins)) {
          lines (x$time, mg[spins[i],], col=rgb(0,i/length(spins),0));
        }
      }
      
    } else {
      if (length (components) == 1) {
        cmp <- rep (components, length (spins));
        
      } else {
        cmp <- components;
      }
      
      if (length(spins) != length(cmp)) {
        stop ("Lengths of spins and components are not compatible!\n", call.=FALSE);

      } else {
        plot (x$time, x$process[(spins[1]-1)*3+cmp[1],], type='l', xlab="Time", ylab="Components", main="Components of the Spins");

        if (length (spins) > 1) {
          for (i in 2:length(spins)) {
            lines (x$time, x$process[(spins[i]-1)*3+cmp[i],], col=rgb(0,i/length(spins),0));
          }
        }
      }
    }
  }

# ----------------------------------------------------------------
component.animation <- function (g=NA, by=10, component=2, time=NA)
  {
    n <- attr (g, "n");  
    T <- length (g$time);
    w <- seq (from=component, by=3, length.out=n);

    if (is.na(time) || time>g$time[T]) {
      t <- T;
    } else {
      t <- length (g$time[g$time<=time]);
    }
    
    tt <- search ();
    if (length(tt["package:rggobi"==tt]) == 0) {
      library (rggobi);
    }
 
    df <- data.frame (x=1:length(w), g$process[w,1]);

    pred <- data.frame(x=1:length(w), y=c(1.3, rep(-1.3, length(w)-1)));
    pre <- ggobi_longitudinal (pred);
    ani <- ggobi_longitudinal (df, g=pre);
    df_a <- ani[1];
 
    for (i in seq (1,t,by=by)) {
      df_a[,2] <- g$process[w,i];
    }
 
    attr (ani, "n") <- n;
    attr (ani, "T") <- t;
    attr (ani, "w") <- w;
    attr (ani, "g") <- g;
    attr (ani, "by") <- by;
    
    return (ani);
  }

replot.animation <- function (ani)
  {
    T <- attr (ani, "T");
    w <- attr (ani, "w");
    by <- attr (ani, "by");
    g <- attr (ani, "g");
    
    df_a <- ani[1];
    for (i in seq (1,T,by=by)) {
      df_a[,2] <- g$process[w,i];
    }
  }

# ----------------------------------------------------------------
rect.demag <- function (a, b, c)
  {
    tmp <- (1/pi)*(((b^2 - c^2)/(2*b*c))*log((sqrt(a^2 + b^2 + c^2) - a)/(sqrt(a^2 + b^2 + c^2) + a)) + ((a^2 - c^2)/(2*a*c))*log((sqrt(a^2 + b^2 + c^2) - b)/(sqrt(a^2 + b^2 + c^2) + b)) + (b/(2*c))*log((sqrt(a^2 + b^2) + a)/(sqrt(a^2 + b^2) - a)) + (a/(2*c))*log((sqrt(a^2 + b^2) + b)/(sqrt(a^2 + b^2) - b)) + (c/(2*a))*log((sqrt(c^2 + b^2) - b)/(sqrt(c^2 + b^2) + b)) + (c/(2*b))*log((sqrt(a^2 + c^2) - a)/(sqrt(a^2 + c^2) + a)) + 2*atan((a*b)/(c*sqrt(a^2 + b^2 + c^2))) + (a^3 + b^3 - 2*c^3)/(3*a*b*c) + (a^2 + b^2 - 2*c^2)*sqrt(a^2 + b^2 + c^2)/(3*a*b*c) + (c/(a*b))*(sqrt(a^2 + c^2) + sqrt(b^2 + c^2)) - ((a^2 + b^2)^(3/2) + (c^2 + b^2)^(3/2) + (a^2 + c^2)^(3/2))/(3*a*b*c));
    tmp
  }

# ---------------------------------------------------
# single spin free energy
spin.energy.config <- function (h.ext=c(0,0,0), K2=0, N.v=c(0,0,0), range=seq(-1.3*pi,1.3*pi,0.02))
  {
    fg1 <- c();

    for (d.theta in range) {
      m <- c(cos(d.theta), sin(d.theta), 0);
      fg1 <- c(fg1, c(d.theta, single.spin.energy(m, N.v, h.ext, K2)));
    }

    fg1 <- t(array (fg1, dim=c(2,length(fg1)/2)));
    plot (fg1, type='l', xlab="Angle", ylab="Energy");
  }
