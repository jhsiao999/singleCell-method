
###  Transformation and Reverse transform functions


## the reverse transform function takes a simplex vector and un-simplexes it on (-infty,infty)

reverse_transform=function(x) 
{
  return(log((abs(x[2:length(x)])+1e-7)/(abs(x[1])+1e-7)));
}

# the transform function simplexes a vector

transform <- function(y) 
{
  temp =c(1,exp(y));
  out=temp/sum(temp);
  return(out)
}
