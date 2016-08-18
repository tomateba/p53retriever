#' Plot the density chart of responsive elements on a DNA sequence
#' 
#'@param p53.table A dataframe containing the responsive elements, as the one returned from the p53.track function. 
#'@param start.at A numeric value indicating the start coordinate of the plot 
#'@param stop.at A numeric value indicating the stop coordinate of the plot
#'@return dp an object of class ggplot
#'@examples
#'data(RE_collection)
#'p53density(RE_collection,-10000,10000)
#'@import ggplot2
#'@export

p53density<-function(p53.table,start.at,stop.at){
  
  bs<-16
  dp<-ggplot(p53.table,aes(start,fill=as.factor(grade),colour=as.factor(grade))) +
    geom_density(alpha=.2) +
    scale_colour_manual(name="",values=c("darkgrey","dodgerblue","green3","indianred","gold2")) +
    scale_fill_manual(name="",values=c("darkgrey","dodgerblue","green3","indianred","gold2")) +
    annotate("point",x = p53.table$start[which(p53.table$grade=="5")], y = 0, alpha = .8, colour="gold2") +
    theme_bw(base_size = bs) +
    scale_x_continuous(limits = c(start.at,stop.at)) +
    labs(y="Density",x="Coordinate")
  return(dp)

}
