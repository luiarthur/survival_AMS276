plotsurv = begin
  R"""
  plotsurv_pch <- function(grid,S,col_l=rep('black',ncol(S)),add=FALSE,addlines=FALSE,...) {
    stopifnot(nrow(S) == length(grid))
    stopifnot(min(S) >= 0 && max(S) <= 1)
    if (!add) plot(grid,xlim=c(-.1,max(grid)),ylim=c(0,1),type='n',...)
    P <- NCOL(S)
    J <- NROW(S)
    for (p in 1:P) {
      if (addlines) lines(grid,S[,p],lwd=1,col=col_l[p],type='s')
      segments(x0=grid[-J], x1=grid[-1], y0=S[-J,p], col=col_l[p],...)
      points(c(grid[-J],grid[-1]), rep(S[-J,p],2), col=col_l[p], pch=c(1,20),...)
      #for (j in 2:J) {
      #  segments(x0=grid[j-1], x1=grid[j], y0=S[j-1,p], col=col_l[p],...)
      #  points(c(grid[j-1],grid[j]), rep(S[j-1,p],2), col=col_l[p], pch=c(1,20),...)
      #}
    }
  }
  """
  R"plotsurv_pch"
end
