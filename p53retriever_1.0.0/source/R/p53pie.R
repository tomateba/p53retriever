#' Plot a pie chart of the grades assigned to responsive elements on a DNA sequence
#' 
#'@param p53.table A dataframe containing the responsive elements, as the one returned from the p53.track function. 
#'@param start.at A numeric value indicating the start coordinate of the site in order to be considered
#'@param stop.at A numeric value indicating the stop coordinate start coordinate of the site in order to be considered
#'@examples
#'data(REs.collection)
#'p53pie(REs.collection,-10000,10000)
#'@export

p53pie<-function(p53.table,start.at,stop.at){

  par(mar=c(5,5,1,1))
  
  p53.table<-p53.table[which((p53.table[,"start"]>start.at)&(p53.table[,"stop"]<stop.at)),]
  slices<-summary(as.factor(p53.table[,"grade"]))[c(1,5,2,4,3)]
  lbls<-c("unlikely","poor","slight","moderate","high")[c(1,5,2,4,3)]
  pct <- round((slices/sum(slices)*100),2)
  lbls <- paste(lbls, pct) # add percents to labels 
  lbls <- paste(lbls,"%",sep="") # ad % to labels
  cols=c("darkgrey","dodgerblue","green3","indianred","gold2")[c(1,5,2,4,3)]
  
pie(slices,col=cols,labels=lbls,radius=0.75,border="white",cex=1.3)
  
}
