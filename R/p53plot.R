#' Plot the position p53 responsive elements on a DNA sequence (the mark is placed at the beginning of the responsive element)
#' 
#'@param p53.table A dataframe containing the responsive elements, as the one returned from the p53.track function. 
#'@param start.at A numeric value indicating the start coordinate of the plot 
#'@param stop.at A numeric value indicating the stop coordinate of the plot
#'@return bp an object of class ggplot
#'@examples
#' data(CDKN1A)
#' hits<-p53track(CDKN1A,seqname="CDKN1A",plot=FALSE)
#' p53plot(hits,0,20000)
#'@import ggplot2
#'@export

p53plot<-function(p53.table,start.at,stop.at){
  bs<-16
  bp<-ggplot(p53.table,aes(start,grade,fill=as.factor(grade))) +
    geom_bar(stat="identity",alpha=.85) +
    scale_fill_manual(name="",values=c("darkgrey","dodgerblue","green3","indianred","gold2")) +
    #annotate("rect",xmin = hits$start, xmax = hits$stop, ymin = 0, ymax = hits$grade, alpha = .8, fill=as.factor(hits$grade)) +
    theme_bw(base_size = bs) +
    scale_x_continuous(limits = c(start.at,stop.at)) +
    scale_y_continuous(limits = c(0,5)) +
    labs(y="Grade",x="Coordinate",title=unique(p53.table$ID))
  return(bp)
}
