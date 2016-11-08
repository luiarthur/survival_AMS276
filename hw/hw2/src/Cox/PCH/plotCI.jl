plotCI = begin
  R"""
  library(rcommon)
  plotCI <- function(grid,S,col_area='cornflowerblue',add=FALSE,...) {
    stopifnot(nrow(S) == length(grid))
    stopifnot(min(S) >= 0 && max(S) <= 1)
    stopifnot(ncol(S) == 2)
    if (!add) plot(grid,xlim=c(-.1,max(grid)),ylim=c(0,1),type='n',fg='grey',...)
    color.btwn(grid, S[,1], S[,2], from=min(grid), to=max(grid), 
               col.area=col_area)
  }
  """
  R"plotCI"
end
