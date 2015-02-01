#' Plot the position p53 responsive elements on a DNA sequence (the mark is placed at the beginning of the responsive element)
#' 
#'@param p53.table A dataframe containing the responsive elements, as the one returned from the p53.track function. 
#'@param start.at A numeric value indicating the start coordinate of the plot 
#'@param stop.at A numeric value indicating the stop coordinate of the plot
#'@export

p53plot<-function(p53.table,start.at,stop.at){

plot(p53.table[which(p53.table[,"grade"]=="1"),"start"],p53.table[which(p53.table[,"grade"]=="1"),"grade"],type="h",lwd=7,ylim=c(0,5),ylab="Grade",xlab="sequence coordinate",main=unique(p53.table[,"ID"]),xlim=c(start.at,stop.at),cex.axis=1.5,cex.lab=1.5,col="darkgrey",frame.plot=F)

points(p53.table[which(p53.table[,"grade"]=="2"),"start"],p53.table[which(p53.table[,"grade"]=="2"),"grade"],type="h",lwd=7,col="dodgerblue")
points(p53.table[which(p53.table[,"grade"]=="3"),"start"],p53.table[which(p53.table[,"grade"]=="3"),"grade"],type="h",lwd=7,col="green3")
points(p53.table[which(p53.table[,"grade"]=="4"),"start"],p53.table[which(p53.table[,"grade"]=="4"),"grade"],type="h",lwd=7,col="indianred")
points(p53.table[which(p53.table[,"grade"]=="5"),"start"],p53.table[which(p53.table[,"grade"]=="5"),"grade"],type="h",lwd=7,col="gold2")

}