#' Plot the density chart of responsive elements on a DNA sequence
#' 
#'@param p53.table A dataframe containing the responsive elements, as the one returned from the p53.track function. 
#'@param start.at A numeric value indicating the start coordinate of the plot 
#'@param stop.at A numeric value indicating the stop coordinate of the plot
#'@examples
#'data(REs.collection)
#'p53density(REs.collection,-10000,10000)
#'@export

p53density<-function(p53.table,start.at,stop.at){

  p53.table<-p53.table[which((p53.table[,"start"]>start.at)&(p53.table[,"stop"]<stop.at)),]
  
  par(mar=c(5,5,1,1))
  
  one<-density(rowMeans(p53.table[which(p53.table[,"grade"]>0),c("start","stop")]))
  two<-density(rowMeans(p53.table[which(p53.table[,"grade"]>1),c("start","stop")]))
  three<-density(rowMeans(p53.table[which(p53.table[,"grade"]>2),c("start","stop")]))
  four<-density(rowMeans(p53.table[which(p53.table[,"grade"]>3),c("start","stop")]))
  five<-density(rowMeans(p53.table[which(p53.table[,"grade"]>4),c("start","stop")]))
  
  ymin<-min(c(one$y,two$y,three$y,four$y,five$y))
  ymax<-max(c(one$y,two$y,three$y,four$y,five$y))
  
  hist(rowMeans(p53.table[which(p53.table[,"grade"]>0),c("start","stop")]),breaks=((stop.at-start.at)/100),xlab="TSS relative coordinate",ylab="Probability density",main=NULL,col="grey85",border="grey85",prob=TRUE,xlim=c(start.at,stop.at), ylim=c((ymin-0.1*(ymax-ymin)),1.2*ymax),cex.axis=1.7,cex.lab=1.7)
  
  #abline(v=0,col="darkgrey",lty=2,lwd=2)
  
  lines(one,col="grey40",lwd=3)
  lines(two,col="dodgerblue",lwd=3)
  lines(three,col="green3",lwd=3)
  lines(four,col="indianred",lwd=3)
  lines(five,col="gold2",lwd=3)
  
  points(p53.table[which(p53.table[,"grade"]=="5"),"start"],rep((ymin-0.05*(ymax-ymin)),length(which(p53.table[,"grade"]=="5"))),col="gold2",cex=1,pch=16)
  
  points(p53.table[which(p53.table[,"grade"]=="4"),"start"],rep((ymin-0.1*(ymax-ymin)),length(which(p53.table[,"grade"]=="4"))),col="indianred",cex=1,pch=16)
  
  legend("top", c(">= 1",">= 2",">= 3",">= 4",">= 5"), fill=c("grey40","dodgerblue","green3","indianred","gold2"), bty="n",cex=1.5,horiz=T,title="Grade")
   
}