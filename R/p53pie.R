#' Plot a pie chart of the grades assigned to responsive elements on a DNA sequence
#' 
#'@param p53.table A dataframe containing the responsive elements, as the one returned from the p53track function. 
#'@param start.at A numeric value indicating the start coordinate of the site in order to be considered
#'@param stop.at A numeric value indicating the stop coordinate start coordinate of the site in order to be considered
#'@return an object of class \code{ggplot}
#'@examples
#'data(RE_collection)
#'p53pie(RE_collection,-10000,10000)
#'@import ggplot2
#'@export

p53pie<-function(p53.table,
                 start.at,
                 stop.at){

  p53.table<-p53.table[which(p53.table$start>start.at & p53.table$stop<stop.at),]
  p53.table$grade<-factor(p53.table$grade,levels=c("1","2","3","4","5"))
  dat <- data.frame(count=as.numeric(table(p53.table$grade)), type=c("unlikely","poor","slight","moderate","high"))
  dat$fraction = dat$count / sum(dat$count)
  dat$ymax = cumsum(dat$fraction)
  dat$ymin = c(0, head(dat$ymax, n=-1))
  dat$label<-paste(dat$count," ",dat$type," (",round(dat$fraction*100,1),"%)",sep="")
  dat$label <- factor(dat$label, levels=as.character(dat$label))
  color.palette <- c("darkgrey","dodgerblue","green3","indianred","gold2")
  
  bs<-16
  pp = ggplot(dat, aes(fill=label, ymax=ymax, ymin=ymin, xmax=4, xmin=0)) +
    geom_rect(colour="white",alpha=0.9) +
    coord_polar(theta="y") +
    xlim(c(0, 4)) +
    #annotate("text",x=0,y=0.5,label=title,size=bs*0.3, fontface="bold") +
    theme_bw(base_size=bs) +
    theme(axis.title=element_blank())+
    theme(panel.grid=element_blank()) +
    theme(axis.text=element_blank()) +
    theme(axis.ticks=element_blank()) +
    theme(panel.border=element_blank())+
    #scale_fill_hue(h=c(90,180),l=c(0,100))
    scale_fill_manual(name="Grade",values=color.palette)
  return(pp)
}
