assign_ld_blocks <- function(l, J){

  nblock <- length(l)
  lstart <- cumsum(c(1, l))[-(nblock + 1)]
  lstop <- lstart + l -1
  #Bookkeeping: Figure out how much/how many replicates of supplied LD we need
  ld_size <- sum(l)
  full_reps <- floor(J/ld_size) # Recall l is list of block sizes
  remainder <- J  - full_reps*ld_size

  block_index <- rep(seq(nblock), full_reps)
  index <- c(rep(seq(ld_size), full_reps)) #, seq(remainder))


  if(full_reps == 0){
    r <- c()
  }else{
    r <- rep(1:full_reps, each = nblock)
  }
  last_block_info <- NULL
  if(remainder > 0){
    index <- c(index, seq(remainder))
    if(l[1] >= remainder){
      blocks_rem <- 0
      final_remainder <- remainder
    }else{
      blocks_rem <- max(which(cumsum(l) <= remainder)) # full blocks in last partial repeat
      final_remainder <- remainder-cumsum(l)[blocks_rem] # partial block in last partial repeat
      block_index <- c(block_index, seq(blocks_rem))
      r <- c(r, rep(full_reps + 1, blocks_rem))
    }
    if(final_remainder > 0){
      block_index <- c(block_index, nblock + 1)
      r <- c(r, full_reps + 1)
      last_block_info <- c(blocks_rem + 1, final_remainder)
      lnew <- c(l, final_remainder)[block_index]
    }else{
      lnew <- l[block_index]
    }
  }else{
    lnew <- l[block_index]
  }
  return(list(block_index = block_index, rep = r,l = lnew, index = index,
              last_block_info = last_block_info))
}
