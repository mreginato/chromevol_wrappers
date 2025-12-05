summaryChromevol <- function(tree, events) {
  require(phangorn)
  tree$tip.label -> spp
  matrix(ncol=2, nrow=length(spp)) -> spp.h
  spp.h[,1] <- spp
  colnames(spp.h) <- c("species", "events")
  
  events -> expecs
  unique(expecs) -> e0
  matrix(ncol=length(e0), nrow=length(spp)) -> spp.m
  rownames(spp.m) <- spp
  colnames(spp.m) <- e0
  spp.m[] <- ""
  
  e0 -> all.events
  sort(all.events) -> all.events
  all.events -> names(all.events)
  all.events
  expecs -> expecs.s
  match(expecs, all.events) -> x
  expecs.s[] <- names(all.events)[x]
  unique(names(all.events)) -> e0
  e0[-match("None",e0)] -> e0
  matrix(ncol=length(e0), nrow=length(spp)) -> spp.m.s
  colnames(spp.m.s) <- e0
  rownames(spp.m.s) <- spp
  spp.m.s[] <- ""
  
  for (i in 1:Ntip(tree)) {
    Ancestors(tree, node=i, type="all") -> n0
    c(sort(n0),i) -> n0
    expecs[n0] -> e0
    e0[which(e0 != "None")] -> e0
    paste(e0, collapse=" | ") -> spp.h[i,2]
    table(e0) -> e0
    match(names(e0), colnames(spp.m)) -> x
    spp.m[i,x] <- e0
    
    expecs.s[n0] -> e0
    e0[which(e0 != "None")] -> e0
    table(e0) -> e0
    match(names(e0), colnames(spp.m.s)) -> x
    spp.m.s[i,x] <- e0
  }
  spp.m.s[spp.m.s[]== ""] <- 0
  mode(spp.m.s) <- "numeric"
  list(history=spp.h, mat=spp.m.s) -> out
  return(out)
}