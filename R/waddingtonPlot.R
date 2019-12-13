#' Plot Waddington Epigenetic Landscape using ggridges
#'
#' @param branches a vector of integer numbers of steady states at each level of the hirarchy tree or developmental path, ordered from top to bottom.
#' @param horizontal_distribution NULL or a list, whose length must be equal with that of vector `branches`, defining relative proportions of widths of the states at the same level. (default: NULL)
#' @param horizontal_distortion NULL or a list of integer vector of length 2, `horizontal_distribution` define proportion of states at the same level, while `horizontal_distortion` define asymmetric curves for each state. (default: NULL)
#' @param vertical_distribution NULL or a integer vecotr, whose length must be 1 less than that of vector `branches`, defining vertical stretch between levels. (default: NULL)
#' @param vertical_distortion NULL or a list of integer vector of length 2. `vertical_distribution` define distances between levlels, while `vertical_distortion` define asymmetric distortion for each betweeness. (default: NULL)
#' @param bottom_width width of bottom of waddington plot
#' @param top_width width of top of waddington plot
#' @param ridge.height height of ridge
#' @param ridge.colors colors of redges
#' @param line.color color of line
#' @param line.type type of line
#' @param line.size size of line
#' @param line.alpha alpha of line
#' @param theme.void logical, whether to use the ggplot2 void theme
#' @param hide.legend logical, whether to hide legend
#' @param do.return logical, whether to return a ggplot object
#'
#' @return NULL if do.return == TRUE. otherwise a ggplot object
#' @export
#'
#' @examples
#' waddingtonPlot(branches = c(1,2,3,4),
#'                horizontal_distribution = list(c(1), c(1,1), c(0.5,0.5,1), c(0.25,0.25,0.5,1)))
waddingtonPlot <- function(branches = c(1,2,4),
                           horizontal_distribution = NULL,
                           horizontal_distortion = NULL,
                           vertical_distribution = NULL,
                           vertical_distortion = NULL,
                           bottom_width = 1000,
                           top_width = 0.3 * bottom_width,
                           ridge.count = 100,
                           ridge.height = 5,
                           ridge.colors = c("red","orange","yellow","green","blue","cyan","purple"),
                           line.color = "black", line.type = 1, line.size = 1, line.alpha = 0.2,
                           theme.void = T, hide.legend = T, do.return = F){
  require("ggplot2")
  require("ggridges")

  # check parameters
  ## horizontal_distribution
  if(is.null(horizontal_distribution)){
    horizontal_distribution <- lapply(branches, function(b) {rep(1, times = b)}) # balanced state
  }else if(is.list(horizontal_distribution) &
           length(horizontal_distribution) == length(branches) &
           all(sapply(horizontal_distribution, length) == branches)){
    pass()
  }else{
    stop("invalid horizontal_distribution setting")
  }
  message('horizontal distribution is ready')

  if(is.null(horizontal_distortion)) {
    horizontal_distortion <- lapply(horizontal_distribution, function(b) {rep(b, each = 2)})
  }else if(length(horizontal_distortion) != 2*sum(branches)){
    stop("elements of horizontal_distortion is not equal with states included in branches")
  }else {
    horizontal_distortion <- unlist(horizontal_distortion)
    cosine_left <- horizontal_distortion[seq(1, length(horizontal_distortion), 2)]
    cosine_right <- horizontal_distortion[seq(2, length(horizontal_distortion), 2)]
    cosine_left <- cosine_left/(cosine_left + cosine_right)
    cosine_right <- cosine_right/(cosine_left + cosine_right)
    
    horizontal_distortion <- list()
    cosine_count <- 0
    for(i in 1:length(horizontal_distribution)) {
      distribution_level <- horizontal_distribution[[i]]
      distortion_level <- c()
      for(j in 1:length(distribution_level)) {
        cosine_count <- cosine_count + j
        distortion_level <- c(distortion_level, 
          cosine_left[cosine_count] * distribution_level[j], 
          cosine_right[cosine_count] * distribution_level[j])
      }
      horizontal_distortion[[i]] <- distortion_level
    }
  }
  message('horizontal distortion is ready')

  

  ## vertical_distribution
  if(is.null(vertical_distribution)){
    vertical_distribution <- rep(1, times = length(branches) - 1)
  }else if(is.vector(vertical_distribution) & length(vertical_distribution) == (length(branches) -1)) {
    pass()
  }else {
    stop("vertical_distribution must be list and its length must be length(branches)-1")
  }
  vertical_distribution <- round(ridge.count * vertical_distribution/sum(vertical_distribution)/2)
  message('vertical distribution is ready')

  if(is.null(vertical_distortion)) {
    vertical_distortion <- lapply(vertical_distribution, function(b) {rep(b, each = 2)})
  }else if(length(unlist(vertical_distortion)) != 2*(length(branches) - 1)){
    stop("elements of vertical_distortion is not 1 less than levels included in branches")
  }else {
    vertical_distortion <- unlist(vertical_distortion)
    cosine_up <- vertical_distortion[seq(1, length(vertical_distortion), 2)]
    cosine_dn <- vertical_distortion[seq(2, length(vertical_distortion), 2)]
    cosine_up <- cosine_up/(cosine_up + cosine_dn)
    cosine_dn <- cosine_dn/(cosine_up + cosine_dn)
    
    vertical_distortion <- list()
    cosine_count <- 0
    for(i in 1:length(vertical_distribution)) {
      distribution_level <- vertical_distribution[[i]]
      distortion_level <- c()
      for(j in 1:length(distribution_level)) {
        cosine_count <- cosine_count + j
        distortion_level <- c(distortion_level, 
          cosine_up[cosine_count] * distribution_level[j], 
          cosine_dn[cosine_count] * distribution_level[j])
      }
      vertical_distortion[[i]] <- round(distortion_level)
    }
  }
  message('vertical distortion is ready')

  
  # skeleton
  skeleton.curves <- lapply(1:length(branches), function(x){
    bias <- horizontal_distribution[[x]]
    sections <- unname(quantile(1:bottom_width, probs = bias/sum(bias), type = 1))
    curves <- sapply(1:(length(sections)/2), function(x) {
      cos.left <- cos((1:sections[2*x-1])*pi/sections[2*x-1])
      cos.right <- cos(pi+ (1:sections[2*x])*pi/sections[2*x])
      return(c(cos.left, cos.right))
    })
    return(unlist(curves)[1:bottom_width])
  }) # get the main wave curves of each one in branches

  # smooth
  ggData <- NULL
  waves_sum <- sum(unlist(vertical_distortion))
  waves_cur <- waves_sum
  ratio <- 1-top_width/bottom_width
  for(i in 1:(length(branches)-1)){
    curve1 <- skeleton.curves[[i]]
    curve2 <- skeleton.curves[[i+1]]
    curve1.5 <- curve1/2 + curve2/2
    # up part
    distortion_up <- vertical_distortion[[i]][1]
    for(j in 0:(distortion_up-1)){
      alpha <- j/distortion_up
      curve <- (1-alpha)*curve1 + alpha*curve1.5
      x_shift <- ratio*(bottom_width/2) * waves_cur/waves_sum
      ggData <- rbind(ggData, data.frame(x = (1:bottom_width) * (1 - ratio * waves_cur / waves_sum) + x_shift,
                                         y = unname(waves_cur), 
                                         shift = curve))
      waves_cur <- waves_cur - 1
    }
    # dn part
    distortion_dn <- vertical_distortion[[i]][2]
    for(j in 0:(distortion_dn-1)){
      alpha <- j/distortion_dn
      curve <- (1-alpha)*curve1.5 + alpha*curve2
      ggData <- rbind(ggData, data.frame(x = (1:bottom_width) * (1 - ratio * waves_cur / waves_sum) + x_shift,
                                         y = unname(waves_cur), 
                                         shift = curve))
      waves_cur <- waves_cur - 1
    }
  }

  # plot
  p <- ggplot() +
    geom_ridgeline(data = ggData,
                   mapping = aes(x, y, height = shift + 2, group = y, fill = y),
                   scale = ridge.height,
                   color = line.color,
                   linetype = line.type,
                   size = line.size,
                   alpha = line.alpha) +
    scale_fill_gradientn(colours = ridge.colors)
  if(theme.void)  p <- p + theme_void()
  if(hide.legend) p <- p + theme(legend.position = "none")


  if(do.return){
    return(p)
  }else{
    print(p)
  }
}
