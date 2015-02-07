#' Locate candidate p53 responsive elements (full and half) on a DNA sequence (working function to detect orphans)
#' 
#' @param seq A character string containing the sequence. The sequence must be composed exclusively of DNA bases (a,c,t,g)
#' @return A dataframe containing the responsive elements located on the input sequence.

#---------------------------------
# Helper function is.between
#---------------------------------

is.between <- function(x, a, b) {
  x > a & x < b
}

#' @export
#' 
#--------------------------------
# Halves and full matrix
#--------------------------------

p53sDebug<-function(seq.ini){

data(bkg.models)
subject<-DNAString(seq.ini)
couples.matrix<-NULL
occurrences.matrix<-NULL
couples.flag<-0
plus.matches <- matchPWM(pwm.1, subject, min.score="80%")     

if(length(plus.matches)>0){
  sequence<-as.character(plus.matches)
  start.seq<-start(plus.matches)
  stop.seq<-end(plus.matches)
  n.mut<-unlist(lapply(sequence, function(x) ((1000-PWMscoreStartingAt(pwm.1, x, starting.at=1))/100)))
  #ww<-unlist(lapply(sequence, function(x) substr(x,5,6)))
  occurrences.matrix<-cbind(sequence,start.seq,stop.seq,n.mut)
  colnames(occurrences.matrix)<-c("seq","start","stop","n.mut")
  
  for(n in 1:nrow(occurrences.matrix)){
    n.couples.run<-which(
      (is.between(as.numeric(occurrences.matrix[,"start"])-as.numeric(occurrences.matrix[n,"start"]),9,24))&
        (as.numeric(occurrences.matrix[,"n.mut"])+as.numeric(occurrences.matrix[n,"n.mut"])<(3)))
    if(length(n.couples.run)>0){
      couples.flag<-1
      first<-t(matrix(rep(occurrences.matrix[n,],length(n.couples.run)),ncol=length(n.couples.run)))
      second<-matrix(occurrences.matrix[n.couples.run,],nrow=length(n.couples.run))
      couples.matrix.run<-cbind(first, second, as.numeric(second[,2])-(as.numeric(first[,3])+1), 
                                as.numeric(second[,4])+as.numeric(first[,4]))
      couples.matrix<-rbind(couples.matrix,couples.matrix.run)
    }
  }
}

#------------------------------------------------------
### half sites
#------------------------------------------------------
halves.matrix<-NULL
resume.matrix.halves<-NULL
if(length(which(occurrences.matrix[,"n.mut"]=="0"))>0) {

halves.matrix<-matrix(occurrences.matrix[occurrences.matrix[,"n.mut"]=="0",],nrow=length(which(occurrences.matrix[,"n.mut"]=="0"))) 
ww.half<-apply(halves.matrix,1,function(x) substr(x[1],5,6))
halves.matrix[,ncol(halves.matrix)]<-ww.half
colnames(halves.matrix)<-c("sequence","start","stop","WW1")

n.mut.tot<-apply(halves.matrix,1,function(x) ((1000-PWMscoreStartingAt(pwm.6, DNAString(x["sequence"]), starting.at=1))/100))

if(length(which(n.mut.tot==0))>0){
  
halves.matrix<-matrix(halves.matrix[which(n.mut.tot==0),],nrow=length(which(n.mut.tot==0)))
colnames(halves.matrix)<-c("sequence","start","stop","WW1")
WW2<-rep("",nrow(halves.matrix))
spacer<-rep(0,nrow(halves.matrix))
n.mut.tot<-rep(0,nrow(halves.matrix))
mutations<-rep("0000000000",nrow(halves.matrix))

halves.matrix.2<-cbind(as.data.frame(halves.matrix),WW2,spacer,mutations,n.mut.tot)

halves.matrix.2[,"start"]<-as.numeric(halves.matrix[,"start"])
halves.matrix.2[,"stop"]<-as.numeric(halves.matrix[,"stop"])
halves.matrix.2[,"spacer"]<-as.numeric(spacer)
halves.matrix.2[,"n.mut.tot"]<-as.numeric(n.mut.tot)
halves.matrix.2<-transform(halves.matrix.2,
                           sequence=as.character(sequence),
                           mutations=as.character(mutations),
                           WW2=as.character(WW2),
                           WW1=as.character(WW1))

halves.matrix.2<-halves.matrix.2[,c("start","stop","spacer","n.mut.tot","sequence","mutations","WW1","WW2")]

grades<-rep(1,length=nrow(halves.matrix.2))
grades[halves.matrix.2[,"WW1"]=="AT"]<-2
labels<-rep("half",length=nrow(halves.matrix.2))

resume.matrix.halves<-cbind(halves.matrix.2,labels,grades)
}
}

#---------------------------------------------------------------------
# Analysis of full sites
#---------------------------------------------------------------------

resume.matrix.full<-class.matrix<-NULL
if(couples.flag==1){
  colnames(couples.matrix)<-c("seq.1","start.1","stop.1","n.mut.1","seq.2","start.2","stop.2","n.mut.2","spacer","n.mut.tot")
  work.matrix<-couples.matrix
  pattern.matrix<-matrix(0,nrow=nrow(work.matrix),ncol=20)
  ww.matrix<-matrix("",nrow=nrow(work.matrix),ncol=2)
  colnames(pattern.matrix)<-c("R1","R2","R3","C1","W1","W2","G1","Y1","Y2","Y3",
                              "R4","R5","R6","C2","W3","W4","G2","Y4","Y5","Y6")
  colnames(ww.matrix)<-c("WW1","WW2")
  double.pwm<-cbind(pwm.6,pwm.6)
  
  for(p in (1:nrow(work.matrix))){
    pair.seq<-paste(couples.matrix[p,"seq.1"],couples.matrix[p,"seq.2"],sep="")
    line<-vector(mode="numeric",length=20)
    for (l in (1:20)){
      letter<-substr(pair.seq,l,l)
      line[l]<-as.numeric(double.pwm[as.character(letter),l]=="0")
    }
    pattern.matrix[p,]<-line
    ww.matrix[p,1]<-as.character(substr(pair.seq,5,6))
    ww.matrix[p,2]<-as.character(substr(pair.seq,15,16))
    
  }
  
  pattern<-as.data.frame(pattern.matrix)
  ww<-as.data.frame(ww.matrix)
  
  couple.matrix<-cbind(work.matrix,pattern,ww) 
  
  sequence<-apply(couple.matrix,1,function(x) toupper(substr(seq.ini,x["start.1"],x["stop.2"])))
  mutations<-apply(couple.matrix,1,function(x) paste(paste(x[11:20],collapse=""),paste(rep("n",length=x["spacer"]),collapse=""),paste(x[21:30],collapse=""),sep=""))
  
  smart.matrix<-cbind(matrix(work.matrix[,c("start.1","stop.2","spacer","n.mut.tot")],nrow=nrow(couple.matrix)),sequence,mutations,ww,pattern)  
  colnames(smart.matrix)[1:2]<-c("start","stop")
  
  smart.matrix[,"n.mut.tot"]<-rowSums(pattern)
  smart.matrix[,"spacer"]<-as.numeric(couples.matrix[,"spacer"])
  smart.matrix[,"start"]<-as.numeric(couples.matrix[,"start.1"])
  smart.matrix[,"stop"]<-as.numeric(couples.matrix[,"stop.2"])

  
class.matrix<-transform(smart.matrix,
                        sequence=as.character(sequence),
                        mutations=as.character(mutations),
                        WW1=as.character(WW1),
                        WW2=as.character(WW2))
grades<-vector(mode="integer",length=nrow(class.matrix))
labels<-vector(mode="character",length=nrow(class.matrix))

}
return(class.matrix)

}
