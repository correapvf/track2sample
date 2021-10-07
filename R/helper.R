ceiling_dec <- function(x, level=3) round(x + 5*10^(-level-1), level)

divide_sample_helper <- function(distance, Ninds, thr_cor, dir, thr, min_thr) {
  if (dir == "F") {
    sample_id <- cumsum(distance) %/% thr_cor
  } else {
    sample_id <- rev(cumsum(rev(distance)) %/% thr_cor)
  }
  df <- data.table::data.table(distance,Ninds,sample_id)
  df <- df[, lapply(.SD, sum), by=.(sample_id), .SDcols=c('distance','Ninds')]
  df <- df[distance >= thr*min_thr]
  ninds <- sum(df$Ninds)
  return(ninds)
}


divide_sample <- function(distance, thr, Ninds, max_thr, min_thr, direction){
  total_disntace <- sum(distance)
  nsamples <- total_disntace %/% thr
  
  if (nsamples == 0) { # if distance is lower than thr
    return("00")
  }
  
  # calculate number of samples
  thr_cor_max <- ceiling_dec(min(total_disntace / nsamples, thr*max_thr+.5))
  thr_cor_min <- ceiling_dec(max(total_disntace / (nsamples+1), thr*min_thr+.5))
  Ninds[is.na(Ninds)] <- 0
  
  opcoes <- expand.grid(thr_cor = c(thr_cor_max, thr_cor_min), dir = direction, ninds = 0)
  
  for (i in 1:nrow(opcoes)) {
    opcoes$ninds[i] <- divide_sample_helper(distance, Ninds, opcoes$thr_cor[i], opcoes$dir[i], thr, min_thr)
  }
  
  # check which orientation will get more abundance
  i <- which.max(opcoes$ninds)
  if (opcoes$dir[i] == "F") {
    sample_id <- cumsum(distance) %/% opcoes$thr_cor[i]
  } else {
    sample_id <- rev(cumsum(rev(distance)) %/% opcoes$thr_cor[i])
  }
  
  
  # return sample id
  return(sprintf("%02d", sample_id))
  
}


create_sample <- function(track, annot, thr, min_ninds, direction, min_thr, max_thr){
  
  # get a sample id by each dive and each habitat
  track2 <- data.table::copy(track)
  track2[, sample_id := sprintf("%03d", data.table::rleid(dive, habitat))]
  
  annot2 <- annot[,.(Ninds = sum(Ninds)), by=.(dive, time)]
  track2[annot2, Ninds_tmp := Ninds, on=.(dive, time)]
  
  # divide samples id by the distance
  track2[, sample := paste(sample_id, divide_sample(distance, thr, Ninds_tmp, max_thr=max_thr, min_thr=min_thr, direction=direction), sep='-'), by=.(sample_id)]
  
  dists <- track2[, .(dists = sum(distance)), by=.(sample)]
  
  annot2 <- track2[annot, .(sample, morphotype, Ninds), on=.(dive, time)]
  
  # summary table
  tmp <- annot2[, .(Nspp = length(unique(morphotype)), Ninds = sum(Ninds)), by=.(sample)]
  tmp[dists, dists := dists, on=.(sample)]
  
  tmp[track2, c('dive','habitat') := .(dive, habitat), on=.(sample)]
  tmp[, keep := dists >= thr*min_thr & Ninds >= min_ninds]
  
  # track and annot table
  track2[, c('sample_id', 'Ninds_tmp') := NULL]
  track2[tmp, keep := keep, on = .(sample)]

  annot2 <- data.table::copy(annot)
  annot2[track2, c('sample', 'keep') := .(sample, keep), on=.(dive, time)]
  
  return(list(track = track2, annot = annot2, summary = tmp))
}

utils::globalVariables(c("." ,".SD", "dive", "time", "distance", "Ninds", "habitat", "morphotype",
                         "keep", "Ninds_tmp", "sample_id"))
