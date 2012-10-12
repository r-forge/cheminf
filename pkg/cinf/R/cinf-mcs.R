# Searching for the Maximum Common Substructures (MCS)
# using the clique search in compatibility graph

make_compatibility_graph <- function(
  connTable1,	# connection table for graph 1
  distMatrix1,	# distance matrix for graph 1
  label1,		# vertex labels for graph 1
  connTable2,	# connection table for graph 2
  distMatrix2,	# distance matrix for graph 2
  label2		# vertex labels for graph 2
) 
{
  size1 <- dim(connTable1)[1]
  size2 <- dim(connTable2)[1]
  index1 <- integer(size1 * size2)
  index2 <- integer(size1 * size2)
  compat_graph_size <- 0
  for (i in 1:size1) {
    for (j in 1:size2) {
	  if (label1[i] == label2[j]) {
	    compat_graph_size <- compat_graph_size + 1
		index1[compat_graph_size] <- i
		index2[compat_graph_size] <- j
	  }
	}
  }
  compatGraph <- matrix(0, compat_graph_size, compat_graph_size)
  index1 <- index1[1:compat_graph_size]
  index2 <- index2[1:compat_graph_size]
  for (i in 1:(compat_graph_size-1)) {
    for (j in (i+1):compat_graph_size) {
	  if ((connTable1[index1[i],index1[j]] == connTable2[index2[i], index2[j]]) &&
	    (distMatrix1[index1[i],index1[j]] == distMatrix2[index2[i], index2[j]]) &&
		(index1[i] != index1[j]) &&
		(index2[i] != index2[j])) {
	    compatGraph[i,j] <- 1
		compatGraph[j,i] <- 1
	  } else {
	    compatGraph[i,j] <- 0
		compatGraph[j,i] <- 0
	  }
	}
  }
  list(compatGraph=compatGraph, index1=index1, index2=index2)
}

# Finds maximum common substructure
find_mcs <- function(
  connTable1,	# connection table for graph 1
  distMatrix1,	# distance matrix for graph 1
  label1,		# vertex labels for graph 1
  connTable2,	# connection table for graph 2
  distMatrix2,	# distance matrix for graph 2
  label2		# vertex labels for graph 2
) 
{

  search_clique <- function(graph) {

    recursive_search_clique <- function(numCandidates, candidates) {
      if (!numCandidates) {
	    # a clique found
	    if (sizeClique > maxSizeClique) {
	      maxSizeClique <<- sizeClique
		  for (i in 1:sizeClique) {
		    maxClique[i] <<- clique[i]
		  }
	    }
	  } else {
	    sizeClique <<- sizeClique + 1
		for (i in 1:numCandidates) {
		  clique[sizeClique] <<- candidates[i]
		  newNumCandidates <- 0
		  newCandidates <- integer(numCandidates - 1)
		  if (i < numCandidates) {
		    for (j in (i+1):numCandidates) {
			  if (graph[candidates[i],candidates[j]]) {
			    newNumCandidates <- newNumCandidates + 1
				newCandidates[newNumCandidates] <- candidates[j]
			  }
			}
		  }
		  if (sizeClique + newNumCandidates > maxSizeClique) {
		    recursive_search_clique(newNumCandidates, newCandidates)
          }		  
		}
	  }
    }

    sizeGraph <- dim(graph)[1]
    clique <- integer(sizeGraph)
    candidates <- integer(sizeGraph)
    maxSizeClique <- 1
    for (i in 1:sizeGraph) {
      clique[1] <- i
	  sizeClique <- 1
	  numCandidates <- 0
	  if (i < sizeGraph) {
	    for (j in (i+1):sizeGraph) {
	      if (graph[i,j]) {
		    numCandidates <- numCandidates + 1
		    candidates[numCandidates] <- j
		    if (sizeClique + numCandidates > maxSizeClique) {
			  cat(sprintf("i=%d j=%d numCandidates=%d\n", i, j, numCandidates)); flush.console()
		      recursive_search_clique(numCandidates, candidates)
		    }
	  	  }
	    }
	  }
    }
  }

  size1 <- dim(connTable1)[1]
  size2 <- dim(connTable2)[1]
  first <- integer(size1)
  second <- integer(size2)
  sizeMCS <- 0
  cat("making compatibility graph...\n"); flush.console()
  res1 <- make_compatibility_graph(connTable1, distMatrix1, label1, connTable2, distMatrix2, label2)
  compatGraph <- res1$compatGraph
  compatGraphSize <- dim(compatGraph)[1]
  cat(sprintf("compatGraphSize=%d\n", compatGraphSize)); flush.console()
  index1 <- res1$index1
  index2 <- res1$index2
  maxClique <- integer(compatGraphSize)
  cat("searching cliques...\n"); flush.console()
  search_clique(compatGraph)
  if (sizeMCS > 1) {
    for (i in 1:sizeMCS) {
	  first[i] <- index1[maxClique[i]]
	  second[i] <- index2[maxClique[i]]
	}
  }  else {
	for (i in 1:sizeMCS) {
	  first[i] <- index1[1]
      second[i]<- index2[1]
    }
  }
}

test_mcs_1 <- function() {
  source("mol2.R")
  cat("reading mdb...\n"); flush.console()
  mdb <- read_mol2("ligands.mol2")
  mol1 <- mdb[[1]]
  mol2 <- mdb[[2]]
  cat("making connection tables...\n"); flush.console()
  mol1_ct <- mol_get_ct(mol1)
  mol2_ct <- mol_get_ct(mol2)
  cat("finding distance matrices...\n"); flush.console() 
  mol1_dm <- calc_distance_matrix(mol1_ct)
  mol2_dm <- calc_distance_matrix(mol2_ct)
  mol1_lab <- mol_get_chelabs(mol1)
  mol2_lab <- mol_get_chelabs(mol2)
  cat("finding mcs...\n"); flush.console()
  res <- find_mcs(mol1_ct, mol1_dm, mol1_lab, mol2_ct, mol2_dm, mol2_lab)
}

#test_mcs_1()

