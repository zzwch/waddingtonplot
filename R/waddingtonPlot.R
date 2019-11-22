#' Plot Waddington Epigenetic Landscape using ggridges
#'
#' @param branches a vector of integer numbers of steady states at each level of the hirarchy tree or developmental path, ordered from top to bottom.
#' @param horizontal_distribution NULL or a list, whose length must be the same with vector `branches`, defining proportions in width of the states on the same level. (default: NULL)
#' @param vertical_distribution NULL or a vector, whose length must be 1 less than vector `branches`, defining vertical stretch between levels. (default: NULL)
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
                           vertical_distribution = NULL,
                           bottom_width = 1000,
                           top_width = 0.3 * bottom_width,
                           ridge.height = 5,
                           ridge.colors = c("red","orange","yellow","green","blue","cyan","purple"),
                           line.color = "black", line.type = 1, line.size = 1, line.alpha = 0.2,
                           theme.void = T, hide.legend = T, do.return = F){
  require("ggplot2")
  require("ggridges")

  # check parameters
  ## horizontal_distribution
  if(is.null(horizontal_distribution)){
    horizontal_distribution <- lapply(branches, function(b) {rep(c(1, 1), times = b)})
  }else if(is.list(horizontal_distribution) &
           length(horizontal_distribution) == length(branches) &
           all(sapply(horizontal_distribution, length) == branches)){
    horizontal_distribution <- lapply(horizontal_distribution, rep, each = 2)
  }else{
    stop("invalid horizontal_distribution setting")
  }

  ## vertical_distribution
  if(is.null(vertical_distribution)){
    vertical_distribution <- rep(20, times = length(branches) - 1)
  }else if(length(vertical_distribution) != (length(branches) -1)){
    stop("length of vertical_distribution must be length(branches)-1")
  }


  waveform.curves <- lapply(1:length(branches), function(x){
    bias <- horizontal_distribution[[x]]
    sections <- unname(quantile(1:bottom_width, probs = bias/sum(bias), type = 1))
    curves <- sapply(1:(length(sections)/2), function(x) {
      cos.left <- cos((1:sections[2*x-1])*pi/sections[2*x-1])
      cos.right <- cos(pi+ (1:sections[2*x])*pi/sections[2*x])
      return(c(cos.left, cos.right))
    })
    return(unlist(curves)[1:bottom_width])
  }) # get the main wave curves of each one in branches

  ggData <- NULL
  waves_sum <- sum(vertical_distribution)
  ratio <- 1-top_width/bottom_width
  for(i in 1:(length(branches)-1)){
    curve1 <- waveform.curves[[i]]
    curve2 <- waveform.curves[[i+1]]
    for(j in 0:(vertical_distribution[i]-1)){
      alpha <- j/vertical_distribution[i]
      curve <- (1-alpha)*curve1 + alpha*curve2 + 1
      waves_cur <- waves_sum - j - sign(i-1)*sum(vertical_distribution[1:(i-1)])
      ggData <- rbind(ggData, data.frame(x = (1:bottom_width)*(1-ratio*waves_cur/waves_sum) + ratio*bottom_width*waves_cur/2/waves_sum,
                                         y = unname(waves_cur), shift = curve))
    }
  }
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
