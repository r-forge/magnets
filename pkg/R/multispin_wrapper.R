# MAKEFLAGS="CFLAGS=-O3" R CMD SHLIB explicit_evol.c

# wrapper for multi-spin with noise case

# init and positions should be nx3 arrays

multispin.evol <- function (init=NULL,
                            positions=NULL,
                            volumes=NULL,
                            fixed.ext=NULL,
                            pos.ext=NULL,
                            vol.ext=NULL,
                            coupling.matrix=NULL,
                            T=10000, dt=0.01,
                            extern.field=0,
                            #uniaxial.anisotropy=0,
                            biaxial.anisotropy=0, # effective rescaled value, K2/(mu*Ms^2)
                            demag.tensor=0,
                            given.noise=TRUE,
                            noise.type="gaussian",
                            sd=1, mean=0, # only useful when given.noise==TRUE
                            dw=NULL,
                            record=TRUE,
                            record.step=1,
                            record.list=NULL,
                            alpha=0.1,
                            nu=0.02,
                            error=10^-10,
                            iter=1000,
                            implicit=FALSE)
  {
    if (is.null(init) || ((is.null(positions) || is.null(volumes)) && is.null(coupling.matrix))) {
      stop ("Some argument is missing!\n", call.=FALSE);
    }

    if (!given.noise && !is.null(dw)) {
      stop ("The noise should not be given!\n", call.=FALSE);
    }

    if (!record && !is.null(record.list)) {
      stop ("The record list should not be given!\n", call.=FALSE);
    }
    
    n <- dim(init)[2];

    # external magnets
    if (is.null(fixed.ext)) {
      me <- 0; ne <- 0;
    } else {
      me <- fixed.ext;
      ne <- dim (fixed.ext)[2];
      if (is.null(coupling.matrix) && length (pos.ext[,1]) != ne) {
        stop ("The number of external magnet positions does not match!\n", call.=FALSE);
      }
    }
    
    # get Cp from positions
    if (!is.null(coupling.matrix)) {
      Cp <- coupling.matrix;
    } else {
      if (ne == 0) {
        Cp <- coupling.matrices (positions, volumes);
      } else {
        if (!is.null(vol.ext) && length (vol.ext) != ne) {
          stop ("The number of external magnet volumes does not match!\n", call.=FALSE);
        }
        pos <- rbind (positions, pos.ext);
        vol <- c(volumes, vol.ext);
        Cp <- coupling.matrices (pos, vol);
      }
    }

    # noise
    if (given.noise) {
      noise.method <- 1;
      if (is.null (dw)) {
        dw <- array (rnorm (n*3*T, mean, sd), dim=c(n*3, T));
      } else if (length (dw) != n*3*T) {
        stop ("The size of noise is not compatible with the system!\n", call.=FALSE);
      }
    } else {
      noise.method <- 0;
      dw <- 0;
    }

    # noise type
    switch (noise.type,
            gaussian = nt <- 0
            );

    # uniaxial & biaxial anisotropy
    #if (length(uniaxial.anisotropy)=1) {
    #  K1 <- rep (uniaxial.anisotropy, n);
    #} else if (length(uniaxial.anisotropy) != n) {
    #  stop ("The uniaxial.anisotropy size is not compatible with the system size!\n", call.=FALSE);
    #}

    if (length(biaxial.anisotropy) == 1) {
      K2 <- rep (biaxial.anisotropy, n);
    } else if (length(biaxial.anisotropy) != n) {
      stop ("The biaxial.anisotropy size is not compatible with the system size!\n", call.=FALSE);
    } else {
      K2 <- biaxial.anisotropy;
    }

    # effective uniaxial.anisotropy should be already included in N.v
    if (length(demag.tensor) == 1) {
      N.v <- rep (demag.tensor, 3*n);
    } else if (length(demag.tensor) == 3) {
      N.v <- rep (demag.tensor, n);
    } else if (length(demag.tensor) != 3*n) {
      stop ("The size of demagnetizing tensor is not compatible with the system size!\n", call.=FALSE);
    } else {
      N.v <- demag.tensor;
    }

    # external field
    if (length (extern.field) == 3) {
      h.ext <- rep (extern.field, n);
      h.type <- 0;
    } else if (length (extern.field) == 3*n) {
      h.ext <- extern.field;
      h.type <- 0;
    } else if (length (extern.field) == 3*T) {
      th <- array (extern.field, dim=c(3, T));
      th1 <- c();
      for (i in 1:n) {
        th1 <- rbind (th1, th);
      }
      h.ext <- as.vector (th1);
    } else if (length (extern.field) == 3*n*T) {
      h.ext <- extern.field;
      h.type <- 1;
    } else {
      stop ("The size of external field is not compatible with the system size!\n", call.=FALSE);
    }

    # record, record.list can override record.step
    if (record == TRUE) {
      rc <- 1;
    } else {
      rc <- 0;
    }
    if (is.null (record.list)) {
      rs <- record.step;
      rl <- 0;
      record.list <- seq (0, T-1, rs);
      rl1 <- length (record.list);
    } else {
      rs <- 0;
      rl <- length (record.list);
      rl1 <- rl;
    }

    if (implicit == TRUE) {
      #so.file <- "implicit_evol.so";
      func <- "implicit_evol";
    } else {
      #so.file <- "explicit_evol.so";
      func <- "explicit_evol";
    }
 
    #dyn.load (so.file); 
    proc <- .C (func,
                as.double (init),
                as.integer (n),
                as.double (me),
                as.integer (ne),
                as.double (Cp),
                as.double (h.ext),
                as.integer (h.type),
                as.double (N.v),
                as.double (dw),
                as.double (K2),
                as.double (alpha),
                as.double (nu),
                as.double (dt),
                as.integer (T),
                as.integer (noise.method),
                as.integer (nt),
                as.integer (rc),
                as.integer (rs),
                as.integer (record.list),
                as.integer (rl),
                as.double (error),
                as.integer (iter),
                process=double(3*n*rl1)
                )$process;
    #dyn.unload (so.file);

    rst <- list (process=array (proc, dim=c(3*n,rl1)), time=(record.list+1));
    class (rst) <- "evolution";
    #attr (rst, "time.slices") <- record.list;
    attr (rst, "call") <- "implicit_evol";
    attr (rst, "dt") <- dt;
    attr (rst, "n") <- n;
    attr (rst, "vol") <- volumes;
    attr (rst, "pos") <- positions;
    #if (h.type == 0) {
    #  attr (rst, "h.ext") <- h.ext;
    #} else {
    #  n3 <- 3 * n;
    #  u <- record.list * n3 + 1;
    #  v <- (record.list+1) * n3;
    #  l <- length (record.list);
    #  u1 <- (1:l-1) * n3 + 1;
    #  v1 <- (1:l) * n3;
    #  tmp <- rep (0, l*n3);
    #  for (i in 1:l) {
    #    tmp[u1[i]:v1[i]] <- h.ext[u[i]:v[i]];
    #  }
    #  attr (rst, "h.ext") <- tmp;
    #}
    attr (rst, "coupling") <- Cp;
    attr (rst, "eff.N") <- N.v; # K1 is already included
    attr (rst, "K2") <- K2;
    attr (rst, "ext.n") <- ne;
    attr (rst, "ext.spins") <- fixed.ext;
    attr (rst, "ext.vol") <- vol.ext;
    attr (rst, "ext.pos") <- pos.ext;

    return (rst);
  }

# -------------------------------------------------------------------------------
# compute the coupling matrices
coupling.matrices <- function (r, V) {
  n <- length (r[,1]);

  # distance vector
  d <- array (0, dim=c(n,n,3));
  if (n > 1) {
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        d[i,j,] <- r[i,]-r[j,];
        d[j,i,] <- -d[i,j,];
      }
    }
  } else {
    d[1,1,] <- c(0,0,0);
  }

  Cp <- array (0, dim=c(n,n,3,3)); # coupling matrices
  for (i in 1:n) {
    for (j in 1:n) {
      dij2 <- sum(d[i,j,]^2);
      Cp[i,j,1,1] <- 3*d[i,j,1]^2 - dij2;
      Cp[i,j,1,2] <- 3*d[i,j,1]*d[i,j,2];
      Cp[i,j,1,3] <- 3*d[i,j,3]*d[i,j,1];
      Cp[i,j,2,1] <- Cp[i,j,1,2];
      Cp[i,j,2,2] <- 3*d[i,j,2]^2 - dij2;
      Cp[i,j,2,3] <- 3*d[i,j,3]*d[i,j,2];
      Cp[i,j,3,1] <- Cp[i,j,1,3];
      Cp[i,j,3,2] <- Cp[i,j,2,3];
      Cp[i,j,3,3] <- 3*d[i,j,3]^2 - dij2;
      if (dij2 != 0) Cp[i,j,,] <- Cp[i,j,,] * V[j] / (4*pi*dij2^(2.5));
    }
  }

  Cp1 <- array (0, dim=c(3*n,3*n));
  for (i in 1:n) {
    for (j in 1:n) {
      i1 <- (i-1)*3 + 1;
      i2 <- i*3;
      j1 <- (j-1)*3 + 1;
      j2 <- j*3;
      Cp1[i1:i2,j1:j2] <- Cp[i,j,,];
    }
  }
  
  return (Cp1);
}
      
