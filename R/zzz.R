.First.lib <- function (lib, pkg) 
{
  library.dynam("Geneland",pkg, lib)
  
  cat('ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo',fill=TRUE)
  cat('o        Geneland is loaded                                       o',fill=TRUE)
  cat('o                                                                 o',fill=TRUE)
  cat('o                * Please *                                       o',fill=TRUE)
  cat('o                                                                 o',fill=TRUE)  
  cat('o        Register on                                              o',fill=TRUE)
  cat('o        http://folk.uio.no/gillesg/Geneland/register.php         o',fill=TRUE)
  cat('o                                                                 o',fill=TRUE)  
  cat('o        See manual on                                            o',fill=TRUE)
  cat('o        http://folk.uio.no/gillesg/Geneland/Geneland.html        o',fill=TRUE)  
  cat('o                                                                 o',fill=TRUE)
  cat('o        This is Geneland-3.1.3                                   o',fill=TRUE)
  cat('o                                                                 o',fill=TRUE)
  cat('ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo',fill=TRUE)
}
