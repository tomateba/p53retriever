#' Plot the position p53 responsive elements on a DNA sequence (the mark is placed at the beginning of the responsive element)
#' 
#'@param p53.table A dataframe containing the responsive elements, as the one returned from the p53track function. 
#'@param start.at A numeric value indicating the start coordinate of the plot 
#'@param stop.at A numeric value indicating the stop coordinate of the plot
#'@return an object of class \code{ggplot}
#'@examples
#' data(CDKN1A)
#' hits<-p53track(CDKN1A,seqname="CDKN1A")
#' p53plot(hits,0,20000)
#'@import ggplot2
#'@export

p53plot<-function(p53.table,start.at,stop.at){
  p53.table$y<-p53.table$grade
  p53.table$grade<-factor(p53.table$grade,levels=c(1,2,3,4,5))
  all.col<-c("darkgrey","dodgerblue","green3","indianred","gold2")
  pres.col<-all.col[which(table(p53.table$grade)>0)]
  selection<-which(p53.table$y>3)
  bs<-16
  bp<-ggplot(p53.table,aes(start,y,fill=grade)) +
    geom_bar(stat="identity",alpha=.85) +
    scale_fill_manual(name="",values=pres.col) +
    annotate("text",x = p53.table[selection,"start"], y = p53.table[selection,"y"], label= p53.table[selection,"label"], alpha = .7, size=bs*0.1, colour="grey20",fontface=3,vjust=-0.1) +
    theme_bw(base_size = bs) +
    scale_x_continuous(limits = c(start.at,stop.at)) +
    scale_y_continuous(limits = c(0,5)) +
    labs(y="Grade",x="Coordinate",title=unique(p53.table$ID))
  bp
  return(bp)
}

#bp2<-bp+geom_vline(xintercept=10000,linetype=2,alpha=0.75,size=bs*0.02,colour="grey30")
#save_plot("cdkn1a.pdf", bp2, base_height = 3.2, base_aspect_ratio=3)
