.First.lib <- function (libname, pkgname)
  {
    library.dynam (chname="magnets", package=pkgname, lib.loc=libname);
  }

.Last.lib <- function (libpath)
  {
    library.dynam.unload (chname="magnets", libpath=libpath);
  }
