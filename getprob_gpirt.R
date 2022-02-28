getprobs_gpirt = function(xs, irfs, thresholds){
  C = ncol(thresholds) - 1
  probs = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(probs) = c("order", "xs","p")
  for (i in 1:length(xs)) {
    ps = matrix(0, nrow=ncol(irfs), ncol=C)
    for (c in 1:C){
      for (iter in 1:ncol(irfs)){
        z1 = thresholds[iter, c] - irfs[i,iter]
        z2 = thresholds[iter, c+1] - irfs[i,iter]
        ps[iter, c] = pnorm(z2)-pnorm(z1)
      }
    }
    ps = colMeans(ps)
    for (c in 1:C){
      probs[nrow(probs) + 1,] = c(c, xs[i],  ps[c])
    }
  }
  return(probs)
}

getprobs_2PL = function(xs, betas){
  C = length(betas)
  probs = data.frame(matrix(ncol = 3, nrow = 0))
  slope = betas[[C]]
  colnames(probs) = c("order", "xs","p")
  for (i in 1:length(xs)) {
    ps = rep(0, C+1)
    ps[C+1] = 1
    for (c in 1:(C-1)){
      lp = betas[[c]] - xs[i]*slope
      ps[c+1] = 1 / (1+ exp(-lp))
    }
    for (c in 1:C) {
      probs[nrow(probs) + 1,] = c(c, xs[i],  ps[c+1]-ps[c]) 
    }
  }
  return(probs)
}
