parseChromevol <- function(res.file = "chromEvol.res", anc.tree.file = "MLAncestralReconstruction.tree",
                           anc.probs.file = "ancestorsProbs.txt", expec.file = "expectations.txt", shifts.file = NULL,
                           thresh = 0.5) {
  require(ape)
  ### res
  readLines(res.file) -> res0
  grep("Previous best ", res0) -> x
  if (length(x) > 0) {
    x-1 -> x
    res0[-c(x:length(res0))] -> res0
  }
  grep("Final optimized likelihood is:", res0) -> s0
  res0[s0:length(res0)] -> res0
  res0[1] -> ln0
  res0[length(res0)] -> aic0
  res0[-c(1,2,length(res0))] -> pars0
  strsplit(ln0, ": ")[[1]][2] -> ln0
  strsplit(aic0, "= ")[[1]][2] -> aic0
  strsplit(pars0, " = ") -> pars0
  do.call(cbind, pars0) -> pars0
  pars0[1,] -> labs0
  pars0[2,] -> pars0
  names(pars0) <- labs0
  c(ln=ln0, AICc=aic0, pars0) -> res
  mode(res) <- "numeric"
  
  ### anc & trees
  read.tree(anc.tree.file) -> t0
  t0$node.label -> anc0
  unlist(lapply(strsplit(anc0, "-"), "[", 2)) -> t0$node.label
  t0 -> tree
  unlist(lapply(strsplit(anc0, "-"), "[", 1)) -> node0
  c((Ntip(t0)+1):(Ntip(t0)+Nnode(t0))) -> nodes.r
  names(nodes.r) <- node0
  read.table(anc.probs.file, header = T, sep="\t", row.names = 1) -> anc0
  sub("-", "", rownames(anc0)) -> rownames(anc0)
  anc0[match(names(nodes.r), rownames(anc0)),] -> anc0
  sub("X", "", colnames(anc0)) -> colnames(anc0)
  round(colSums(anc0),3) -> x
  anc0[,x > 0] -> ancs
  
  ### expecs
  
  target <- "#ALL EVENTS EXPECTATIONS PER NODE"
  end = "#"
  names(nodes.r) -> nodes.key
  names(nodes.key) <- nodes.r
  unlist(lapply(strsplit(tree$tip.label, "-"), "[", 1)) -> tips0
  c(1:length(tips0)) -> names(tips0)
  c(tips0,nodes.key) -> nodes.key
  readLines(expec.file) -> all.events
  all.events[grep(target, all.events):length(all.events)] -> all.events
  all.events[-1] -> all.events
  all.events[1] -> head0
  all.events[2:grep(end, all.events)[1]] -> all.events
  all.events[-length(all.events)] -> all.events
  strsplit(all.events, "\t") -> all.events
  do.call(rbind, all.events) -> all.events
  all.events[,1] -> rownames(all.events)
  all.events[,-1] -> all.events
  mode(all.events) <- "numeric"
  all.events[,-ncol(all.events)] -> all.events
  head0
  colnames(all.events) <- c("Gain", "Loss", "Duplication", "Demi-duplication", "Transitions-base-number")
  
  ### include root
  match(nodes.key, rownames(all.events)) -> x
  which(is.na(x)) -> root0
  rep(0, ncol(all.events)) -> x
  names(x) <- colnames(all.events)
  rbind(all.events, x) -> all.events
  rownames(all.events)[nrow(all.events)] <- nodes.key[root0]
  all.events[match(nodes.key, rownames(all.events)),] -> all.events
  
  ### Get state
  
  names(nodes.key) -> all.nodes.events
  all.nodes.events -> names(all.nodes.events)
  all.nodes.events[] <- NA
  
  ### transitions base number
  
  which(all.events[,5] > 0.5) -> x
  all.nodes.events[x] <- "Transitions-base-number"
  
  ### demi & dupli
  
  all.events[,c(3,4)] -> x
  for (i in 1:nrow(x)) {
    x[i,] -> x0
    x0[which.max(x0)] -> x0
    if (x0 > thresh) {
      paste(all.nodes.events[i], names(x0), sep = " | ") -> all.nodes.events[i]
    }
  }
  
  ### gain & loss
  
  all.events[,c(1,2)] -> x
  for (i in 1:nrow(x)) {
    x[i,] -> x0
    x0[which.max(x0)] -> x0
    if (x0 > thresh) {
      paste(all.nodes.events[i], names(x0), sep = " | ") -> all.nodes.events[i]
    }
  }
  
  ### none
  gsub("NA | ", "", all.nodes.events, fixed=T) -> all.nodes.events
  all.nodes.events[which(is.na(all.nodes.events))] <- "None"
  
  ### shifts 
  
  read.tree(shifts.file) -> t0
  t0$node.label -> anc0
  unlist(lapply(strsplit(anc0, "-"), "[", 2)) -> t0$node.label
  t0 -> shifts
  
  ### return
  list(res=res, all.events=all.events, events=all.nodes.events, anc=ancs, tree=tree, shifts=shifts) -> out
  return(out)
  
}
