"rdiscr" <-
function(f)
  {
    F <- cumsum(f)
    indic <- F<runif(1)
    1+sum(indic)
  }

