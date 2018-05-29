.onAttach <- function(libname, pkgname)
   packageStartupMessage("glmmEP 1.0 loaded.\nCopyright M.P. Wand and J.F.C. Yu 2018.\nFor details on the use of glmmEP, issue the command:\nglmmEPvignette()")

.onUnload <- function(libpath)
    library.dynam.unload("glmmEP",libpath)
