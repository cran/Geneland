.First.lib <- function (lib, pkg) 
{
  library.dynam("Geneland",pkg, lib)

  
  cat('*********************************************************************************',fill=T)
  cat('*   Geneland-1.0.1 is loaded                                                    *',fill=T)
  cat('*                                                                               *',fill=T)
  cat('*   Use help(Geneland) for a quick overview.                                    *',fill=T)
  cat('*                                                                               *',fill=T)
  cat('*   Please have a look at the warning section at:                               *',fill=T)
  cat('*   www.inapg.inra.fr/ens_rech/mathinfo/personnel/guillot/Geneland.html         *',fill=T)
  cat('*                                                                               *',fill=T)     
  cat('*********************************************************************************',fill=T)

}


