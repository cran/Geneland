.First.lib <- function (lib, pkg) 
{
  library.dynam("Geneland",pkg, lib)

  
  cat('*********************************************************************************',fill=T)
  cat('*   Geneland-1.0.7 is loaded                                                    *',fill=T)
  cat('*                                                                               *',fill=T)
  cat('*   Use help(Geneland) for a quick overview.                                    *',fill=T)
  cat('*                                                                               *',fill=T)
  cat('*   Please                                                                      *',fill=T)
  cat('*      - have a look at the warning section at:                                 *',fill=T)
  cat('*        www.inapg.inra.fr/ens_rech/mathinfo/personnel/guillot/Geneland.html    *',fill=T)
  cat('*      - send a short email at guillot@inapg.inra.fr                            *',fill=T)
  cat('*        to notify that you are using Geneland                                  *',fill=T)
  cat('*                                                                               *',fill=T)
  cat('*********************************************************************************',fill=T)

}


