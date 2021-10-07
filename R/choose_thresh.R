
#' Choose and apply a distance threshold to convert a track to samples
#' 
#' \code{choose_thresh} will subdivide a ROV track in segments using multiple distance thresholds
#' For each threshold abundance and richness are calculated.
#' \code{track2segment} will return the track or annotations with a unique sample id based on a threshold.
#' 
#' track data.frame should have columns:
#' \itemize{
#'   \item \code{dive} - char or factor, indicating the current ROV track
#'   \item \code{time} - char, numerical or posixt, current time stamp
#'   \item \code{habitat} - char, numerical or factor, indicating the current habitat. Set to a single value to treat all track as a single habitat (basically ignoring this field).
#'   \item \code{distance} - how much meters (or other units) the ROV moved from the last row. The first row should be 0.
#'  }
#' 
#' annot data.frame should have columns:
#' \itemize{
#'   \item \code{dive} - char or factor, indicating the current ROV track
#'   \item \code{time} - char or posixt, current time stamp
#'   \item \code{morphotype} - char or factor, indicating the current habitat. Set to a single value to treat all track as a single habitat (basically ignoring this field).
#'   \item \code{Ninds} - how much meters (or other units) the ROV moved from the last row. The first row should be 0.
#'  }
#' 
#' The \code{min_thr} and \code{max_thr} are multipliers that allow threshold to be more flexible. 
#' E.g. if a threshold = 100, min_thr = 0.8 and max_thr = 1.1, than the threshold for the segments can vary from 80 (100 x 0.8) and 110 (100 x 1.1).
#' Adding such flexibility allows more observations to be included in the analysis (e.g. 95 m long segment will still be considered as valid, 
#' or you can split 612 m long track in 6 segments of 102, instead of using exactly 100 m and letting out those 12 m).
#' When setting min_thr or max_thr to values different than 1, you should test if differences in segment sizes do not effect abundance and richness.
#' 
#' Time stamps from \code{annot} data.frame must be present in the \code{track} data.frame, or else these annotations will be disregarded.
#' \code{dive} and \code{time} columns should be of the same type in both data.frames.
#' @examples 
#' # Example with dummy track and annotations
#' set.seed(12345)
#' track <- data.frame(
#'     dive = rep(c('da','db'), each=1000),
#'     time = rep(1:1000, 2),
#'     habitat = rep(rep(c('ha','hb'), 7), c(20,300,400,80,110,290,50,130,60,200,160,130,40,30)),
#'     distance = runif(2000))
#' annot <- data.frame(
#'     dive = rep(c('da','db'), each=300),
#'     time = sample(1:1000, 600, replace=TRUE),
#'     morphotype = sample(letters[1:10], 600, replace=TRUE),
#'     Ninds = round(rpois(600, 0.5)+runif(600)))
#' 
#' result <- choose_thresh(track, annot, min_ninds = 20, distances = seq(20, 100, 10))
#' annot_samples <- track2segment(track, annot, 50, return_type = "annot")
#' @param track A data.frame with the moving distance in each n-seconds. 
#' @param annot A data.frame with the video annotations.
#' @param min_ninds Segments below this number of individuals will be rejected.
#' @param distances A vector of distance thresholds to use as sample lengths.
#' @param direction One of \code{"F"} - segments are placed at the start of the track, \code{"R"} - segments are placed at the end of the track, or
#' \code{"both"} - either F or R in order to maximize abundance.
#' @param min_thr,max_thr Minimal and maximal multipliers applied to the distances that creates a tolerance in the final segment sizes.
#' @param add_plot logical, Should plot results for abundance?
#' @rdname track2segment
#' @name track2segment
#' @export
choose_thresh <- function(track, annot, min_ninds = 20, distances = seq(20, 300, 10), direction = "F", min_thr = 1.0, max_thr = 1.0, add_plot = TRUE) {
  track <- data.table::as.data.table(track)
  annot <- data.table::as.data.table(annot)
  
  if (class(track$dive) != class(annot$dive)) stop("'dive' column should be the same type in both 'track' and 'annot' data.frames.")
  if (class(track$time) != class(annot$time)) stop("'time' column should be the same type in both 'track' and 'annot' data.frames.")
  if (min_thr > 1) stop("'min_tr' should be below than 1")
  if (max_thr < 1) stop("'max_tr' should be higher than 1")
  if (!(direction %in% c('F','R'))) direction <- c('F','R')
  
  # iterate for each distance threshold
  results <- data.frame(thr = distances, # the thr used
                       min_richness = 0, # min richness of the samples
                       max_richness = 0, # max richness of the samples
                       richness = 0, # overall richness
                       min_inds = 0, # minimal inds in samples
                       max_inds = 0, # maximum inds in samples
                       sum_inds = 0, # number of inds
                       n_samples = 0, # number of samples
                       distance = 0, # distance that was considered in the samples
                       n_habitat = 0) # number of habitat covered by the segments
  
  pb = utils::txtProgressBar(max = length(distances), style = 3)
  
  for (i in seq_along(distances)) { # for each thr
    # get samples ids
    samples <- create_sample(track, annot, distances[i], min_ninds, direction, min_thr, max_thr)
    
    # select samples above thr
    tmp <- samples[[3]][keep == TRUE]
    
    # store results
    results$max_richness[i] <- max(tmp$Nspp)
    results$min_richness[i] <- min(tmp$Nspp)
    results$richness[i] <- length(unique(samples[[2]][keep == TRUE, morphotype]))
    results$n_samples[i] <- nrow(tmp)
    results$sum_inds[i] <- sum(tmp$Ninds)
    results$min_inds[i] <- min(tmp$Ninds)
    results$max_inds[i] <- max(tmp$Ninds)
    results$distance[i] <- sum(tmp$dists)
    results$n_habitat[i] <- sum(table(tmp$habitat) > 0)
    
    utils::setTxtProgressBar(pb, i)
  }
  
  close(pb)
  if (add_plot) plot(results$thr, results$sum_inds, type='b', xlab="Segment Size", ylab="Total abundance")
  return(results)
}



#' @param thr numerical, a single distance threshold to create samples
#' @param return_type character. What should the function return? \code{track} - track data.frame,
#' \code{annot} -  annot data.frame, \code{summary} - a summary data.frame of each sample, or \code{all} - a list with the three data.frames.
#' @importFrom data.table :=
#' @export
track2segment <- function(track, annot, thr, min_ninds = 20, direction = "F", min_thr = 1.0, max_thr = 1.0, return_type = "annot"){
  track <- data.table::as.data.table(track)
  annot <- data.table::as.data.table(annot)
  
  samples <- create_sample(track, annot, thr, min_ninds, direction, min_thr, max_thr)
  
  out <- switch (return_type,
    "track" = samples[[1]],
    "annot" = samples[[2]],
    "summary" = samples[[3]],
    samples
  )
  
  return(out)
}
