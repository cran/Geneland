.First.lib <- function(package,chname)
  {
    if(.Platform$OS == "unix") library.dynam(chname="Geneland.so",package="Geneland") 
    if(.Platform$OS == "windows") library.dynam(chname="Geneland.dll",package="Geneland") 
  }
 
