#' Locate candidate p53 responsive elements (full and half) on a DNA sequence
#' 
#' @param seq.ini A character string containing the sequence. The sequence must be composed exclusively of DNA bases (a,c,t,g)
#' @return A dataframe containing the responsive elements located on the input sequence.
#' @export 

# Halves and full matrix
#--------------------------------

p53sf<-function(seq.ini){

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

if(length(which(n.mut.tot<1))>0){
  
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

resume.matrix.full<-NULL
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

#-----------------------
# Rules for full sites
#-----------------------

#1) O
#----------------------------------------------------
labels[which(class.matrix[,"n.mut.tot"]==0)]<-"O"
#-----------------------------------------------------

grades[which(class.matrix[,"n.mut.tot"]==0)]<-2

grades[which((class.matrix[,"n.mut.tot"]==0)
             &(class.matrix[,"spacer"]%in%c(1,2))
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-3

grades[which((class.matrix[,"n.mut.tot"]==0)
             &(class.matrix[,"spacer"]==0))]<-4

grades[which((class.matrix[,"n.mut.tot"]==0)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-5

#2) A

#----------------------------------------------------
mut.cols<-c("C1","G2","C2","G1")
labels[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1))]<-"A"
#----------------------------------------------------

mut.cols<-c("C1","G2")

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]==0))]<-3

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-4

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]!=0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-2

#------------------------------------

mut.cols<-c("C2","G1")

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]==0))]<-2

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-3


#3) B
#----------------------------------------------------
mut.cols<-c("W1","W2","W3","W4")
labels[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1))]<-"B"
#----------------------------------------------------

mut.cols<-c("W1","W4")

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]==0))]<-3

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-4

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]%in%c(1,2)))]<-2

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]>2)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-2

#------------------------------------

mut.cols<-c("W2","W3")

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]==0))]<-3

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-4

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]>0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-2


#4) C

#----------------------------------------------------
mut.cols<-c("R1","R4","Y3","Y6")
labels[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1))]<-"C"
#----------------------------------------------------

mut.cols<-c("R1","Y6")

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]==0))]<-4

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-5

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]%in%c(1,2)))]<-2

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]%in%c(1,2))
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-3


grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]>2)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-2

#------------------------------------

mut.cols<-c("R4","Y3")

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]==0))]<-4

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-5

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]%in%c(1,2)))]<-2

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]%in%c(1,2))
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-3


grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]>2)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-2

#------------------------------------



#5) D


#----------------------------------------------------
mut.cols<-c("R2","R3","R5","R6","Y1","Y2","Y4","Y5")
labels[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1))]<-"D"
#----------------------------------------------------

mut.cols<-c("R2","R3","Y4","Y5")

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]==0))]<-3

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-4

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]>0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-2

#------------------------------------

mut.cols<-c("R5","R6","Y1","Y2")

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]==0))]<-2

grades[which((class.matrix[,"n.mut.tot"]==1)
             &(rowSums(class.matrix[,mut.cols])==1)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-4

#------------------------------------



#6) AA

#----------------------------------------------------
mut.cols<-c("C1","G1","C2","G2")
labels[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2))]<-"AA"
#----------------------------------------------------
# same half-site

mut.cols<-c("C1","G1","C2","G2")

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &(rowSums(class.matrix[,mut.cols[1:2]])!=1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &(rowSums(class.matrix[,mut.cols[1:2]])!=1)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-2

# different half-site
#----------------------------------------------------

mut.cols<-c("C1","G1","C2","G2")

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &(rowSums(class.matrix[,mut.cols[1:2]])==1))]<-1

#----------------------------------------------------


# 7) BB

#----------------------------------------------------
mut.cols<-c("W1","W2","W3","W4")
labels[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2))]<-"BB"
#----------------------------------------------------
# same half site

mut.cols<-c("W1","W2","W3","W4")

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &(rowSums(class.matrix[,mut.cols[1:2]])!=1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &(rowSums(class.matrix[,mut.cols[1:2]])!=1)
             &(class.matrix[,"spacer"]==0))]<-3

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &(rowSums(class.matrix[,mut.cols[1:2]])!=1)
             &(class.matrix[,"spacer"]>0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-2

#----------------------------------------------------
#different half site

mut.cols<-c("W1","W2","W3","W4")

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &(rowSums(class.matrix[,mut.cols[1:2]])==1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &(rowSums(class.matrix[,mut.cols[1:2]])==1)
             &(class.matrix[,"spacer"]==0))]<-2



#----------------------------------------------------



# 8) AB

#----------------------------------------------------------
mut.cols.1<-c("C1","G1","C2","G2")
mut.cols.2<-c("W1","W2","W3","W4")
labels[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1))]<-"AB"
#-----------------------------------------------------------
# same half site
mut.cols.1<-c("C1","G1","C2","G2")
mut.cols.2<-c("W1","W2","W3","W4")

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:2])])!=1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:2])])!=1)
             &(class.matrix[,"spacer"]==0))]<-3

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:2])])!=1)
             &(class.matrix[,"spacer"]>0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-2

#----------------------------------------------------
#different half site

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:2])])==1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:2])])==1)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]==0))]<-2

#----------------------------------------------------
# 1st or 4th quarter

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1,4)],mut.cols.2[c(1,4)])])==2)))]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1,4)],mut.cols.2[c(1,4)])])==2))
             &(class.matrix[,"spacer"]==0))]<-3

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1,4)],mut.cols.2[c(1,4)])])==2))
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]==0))]<-4

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1,4)],mut.cols.2[c(1,4)])])==2))
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]>0))]<-2

#----------------------------------------------------
# all 1st or all 4th quarter

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1)],mut.cols.2[c(1)])])==2)|(rowSums(class.matrix[,c(mut.cols.1[c(4)],mut.cols.2[c(4)])])==2))
)]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1)],mut.cols.2[c(1)])])==2)|(rowSums(class.matrix[,c(mut.cols.1[c(4)],mut.cols.2[c(4)])])==2))
             &(class.matrix[,"spacer"]==0))]<-3

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1)],mut.cols.2[c(1)])])==2)|(rowSums(class.matrix[,c(mut.cols.1[c(4)],mut.cols.2[c(4)])])==2))
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]==0))]<-4

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1)],mut.cols.2[c(1)])])==2)|(rowSums(class.matrix[,c(mut.cols.1[c(4)],mut.cols.2[c(4)])])==2))
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]%in%c(1,2)))]<-3

#------------------------------------------------------
# 9) BD

#----------------------------------------------------
mut.cols.1<-c("W1","W2","W3","W4")
mut.cols.2<-c("R2","R3","Y1","Y2","R5","R6","Y4","Y5")
labels[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1))]<-"BD"
#----------------------------------------------------

# same half site
mut.cols.1<-c("W1","W2","W3","W4")
mut.cols.2<-c("R2","R3","Y1","Y2","R5","R6","Y4","Y5")

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:4])])!=1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:4])])!=1)
             &(class.matrix[,"spacer"]==0))]<-2

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:4])])!=1)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]==0))]<-3

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:4])])!=1)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]>0))]<-2

#----------------------------------------------------
#different half site

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:4])])==1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:4])])==1)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]==0))]<-2

#----------------------------------------------------
# 1st or 4th quarter

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1,4)],mut.cols.2[c(1,2,7,8)])])==2)))]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1,4)],mut.cols.2[c(1,2,7,8)])])==2))
             &(class.matrix[,"spacer"]==0))]<-2

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1,4)],mut.cols.2[c(1,2,7,8)])])==2))
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]==0))]<-3

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1,4)],mut.cols.2[c(1,2,7,8)])])==2))
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]>0))]<-2

#----------------------------------------------------
# All 1st or all 4th quarter

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1)],mut.cols.2[c(1,2)])])==2)|(rowSums(class.matrix[,c(mut.cols.1[c(4)],mut.cols.2[c(7,8)])])==2))
)]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1)],mut.cols.2[c(1,2)])])==2)|(rowSums(class.matrix[,c(mut.cols.1[c(4)],mut.cols.2[c(7,8)])])==2))
             &(class.matrix[,"spacer"]==0)
)]<-3

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1)],mut.cols.2[c(1,2)])])==2)|(rowSums(class.matrix[,c(mut.cols.1[c(4)],mut.cols.2[c(7,8)])])==2))
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]==0)
)]<-4

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1)],mut.cols.2[c(1,2)])])==2)|(rowSums(class.matrix[,c(mut.cols.1[c(4)],mut.cols.2[c(7,8)])])==2))
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]>0)
)]<-2

#--------------------------------------------------


# 10) AD

#----------------------------------------------------
mut.cols.1<-c("C1","G1","C2","G2")
mut.cols.2<-c("R2","R3","Y1","Y2","R5","R6","Y4","Y5")
labels[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1))]<-"AD"
#----------------------------------------------------

# same half site
mut.cols.1<-c("C1","G1","C2","G2")
mut.cols.2<-c("R2","R3","Y1","Y2","R5","R6","Y4","Y5")

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:4])])!=1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:4])])!=1)
             &(class.matrix[,"spacer"]==0))]<-2

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:4])])!=1)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-3

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:4])])!=1)
             &(class.matrix[,"spacer"]>0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-2

#----------------------------------------------------
#different half site

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:4])])==1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:4])])==1)
             &(class.matrix[,"spacer"]==0))]<-2

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:4])])==1)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]==0))]<-3


#----------------------------------------------------
# 1st or 4th quarter

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1,4)],mut.cols.2[c(1,2,7,8)])])==2)))]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1,4)],mut.cols.2[c(1,2,7,8)])])==2))
             &(class.matrix[,"spacer"]==0))]<-2

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1,4)],mut.cols.2[c(1,2,7,8)])])==2))
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]==0))]<-3

#--------------------------------------------------

# All 1st or all 4th quarter

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1)],mut.cols.2[c(1,2)])])==2)|(rowSums(class.matrix[,c(mut.cols.1[c(4)],mut.cols.2[c(7,8)])])==2))
)]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1)],mut.cols.2[c(1,2)])])==2)|(rowSums(class.matrix[,c(mut.cols.1[c(4)],mut.cols.2[c(7,8)])])==2))
             &(class.matrix[,"spacer"]==0)
)]<-3

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1)],mut.cols.2[c(1,2)])])==2)|(rowSums(class.matrix[,c(mut.cols.1[c(4)],mut.cols.2[c(7,8)])])==2))
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]==0)
)]<-4

#-------------------------------------------------------------------------------------


# 11) DD

#----------------------------------------------------

mut.cols<-c("R2","R3","Y1","Y2","R5","R6","Y4","Y5")
labels[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2))]<-"DD"
#----------------------------------------------------

# same half site
mut.cols<-c("R2","R3","Y1","Y2","R5","R6","Y4","Y5")


grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &(rowSums(class.matrix[,mut.cols[1:4]])!=1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &(rowSums(class.matrix[,mut.cols[1:4]])!=1)
             &(class.matrix[,"spacer"]==0))]<-2

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &(rowSums(class.matrix[,mut.cols[1:4]])!=1)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]==0))]<-3


#----------------------------------------------------
#different half site

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &(rowSums(class.matrix[,mut.cols[1:4]])==1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &(rowSums(class.matrix[,mut.cols[1:4]])==1)
             &(class.matrix[,"spacer"]==0))]<-2

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &(rowSums(class.matrix[,mut.cols[1:4]])==1)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]==0))]<-3

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &(rowSums(class.matrix[,mut.cols[1:4]])==1)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]%in%c(1,2)))]<-2


#----------------------------------------------------
# 1st or 4th quarter

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &(rowSums(class.matrix[,mut.cols[c(1,2,7,8)]])==2))]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &(rowSums(class.matrix[,mut.cols[c(1,2,7,8)]])==2)
             &(class.matrix[,"spacer"]==0))]<-2

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &(rowSums(class.matrix[,mut.cols[c(1,2,7,8)]])==2)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]==0))]<-3


grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &(rowSums(class.matrix[,mut.cols[c(1,2,7,8)]])==2)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]%in%c(1,2)))]<-2

#--------------------------------------------------
# All 1st or all 4th quarter


grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &((rowSums(class.matrix[,mut.cols[c(1,2)]])==2)|(rowSums(class.matrix[,mut.cols[c(7,8)]])==2))
)]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &((rowSums(class.matrix[,mut.cols[c(1,2)]])==2)|(rowSums(class.matrix[,mut.cols[c(7,8)]])==2))
             &(class.matrix[,"spacer"]==0))]<-3

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &((rowSums(class.matrix[,mut.cols[c(1,2)]])==2)|(rowSums(class.matrix[,mut.cols[c(7,8)]])==2))
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]==0))]<-4


#--------------------------------------------------
# 12) CC
#----------------------------------------------------

mut.cols<-c("R1","Y3","R4","Y6")
labels[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2))]<-"CC"
#----------------------------------------------------

mut.cols<-c("R1","Y3","R4","Y6")

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2))]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &(class.matrix[,"spacer"]==0))]<-4

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]==0))]<-5

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &(class.matrix[,"spacer"]%in%c(1,2)))]<-2

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols])==2)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]%in%c(1,2)))]<-3

#--------------------------------------------------


# 13) BC

#----------------------------------------------------
mut.cols.1<-c("W1","W2","W3","W4")
mut.cols.2<-c("R1","Y3","R4","Y6")
labels[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1))]<-"BC"

#----------------------------------------------------
# same half site

mut.cols.1<-c("W1","W2","W3","W4")
mut.cols.2<-c("R1","Y3","R4","Y6")

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:2])])!=1)
)]<-1


grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:2])])!=1)
             &(class.matrix[,"spacer"]==0))]<-2

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:2])])!=1)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]==0))]<-3


#----------------------------------------------------
#different half site

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:2])])==1)
)]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:2])])==1)
             &(class.matrix[,"spacer"]==0))]<-3

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:2])])==1)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]==0))]<-4

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:2])])==1)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]%in%c(1,2)))]<-2

#----------------------------------------------------
# all 1st or all 4th quarter

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1)],mut.cols.2[c(1)])])==2)|(rowSums(class.matrix[,c(mut.cols.1[c(4)],mut.cols.2[c(4)])])==2))
)]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1)],mut.cols.2[c(1)])])==2)|(rowSums(class.matrix[,c(mut.cols.1[c(4)],mut.cols.2[c(4)])])==2))
             &(class.matrix[,"spacer"]==0))]<-3

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1)],mut.cols.2[c(1)])])==2)|(rowSums(class.matrix[,c(mut.cols.1[c(4)],mut.cols.2[c(4)])])==2))
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]==0))]<-4

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1)],mut.cols.2[c(1)])])==2)|(rowSums(class.matrix[,c(mut.cols.1[c(4)],mut.cols.2[c(4)])])==2))
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]%in%c(1,2)))]<-2

#--------------------------------------------------

# 14) CD

#----------------------------------------------------
mut.cols.1<-c("R1","Y3","R4","Y6")
mut.cols.2<-c("R2","R3","Y1","Y2","R5","R6","Y4","Y5")
labels[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1))]<-"CD"
#----------------------------------------------------

mut.cols.1<-c("R1","Y3","R4","Y6")
mut.cols.2<-c("R2","R3","Y1","Y2","R5","R6","Y4","Y5")

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
)]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(class.matrix[,"spacer"]==0))]<-3

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]==0))]<-4


grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]%in%c(1,2)))]<-2

#----------------------------------------------------

# 15) AC
#----------------------------------------------------
mut.cols.1<-c("C1","G1","C2","G2")
mut.cols.2<-c("R1","Y3","R4","Y6")
labels[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1))]<-"AC"
#----------------------------------------------------

# same half site
mut.cols.1<-c("C1","G1","C2","G2")
mut.cols.2<-c("R1","Y3","R4","Y6")

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:2])])!=1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:2])])!=1)
             &(class.matrix[,"spacer"]==0))]<-2

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:2])])!=1)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-3

#----------------------------------------------------
#different half site


grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:2])])==1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:2])])==1)
             &(class.matrix[,"spacer"]==0))]<-2

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,c(mut.cols.1[1:2],mut.cols.2[1:2])])==1)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-3

#----------------------------------------------------
# all 1st or all 4th quarter

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1)],mut.cols.2[c(1)])])==2)|(rowSums(class.matrix[,c(mut.cols.1[c(4)],mut.cols.2[c(4)])])==2))
)]<-1

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1)],mut.cols.2[c(1)])])==2)|(rowSums(class.matrix[,c(mut.cols.1[c(4)],mut.cols.2[c(4)])])==2))
             &(class.matrix[,"spacer"]==0))]<-3

grades[which((class.matrix[,"n.mut.tot"]==2)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &((rowSums(class.matrix[,c(mut.cols.1[c(1)],mut.cols.2[c(1)])])==2)|(rowSums(class.matrix[,c(mut.cols.1[c(4)],mut.cols.2[c(4)])])==2))
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-4

#--------------------------------------------------

### 16) A+C+C

#----------------------------------------------------
mut.cols.1<-c("C1","G1","C2","G2")
mut.cols.2<-c("R1","Y3","R4","Y6")
labels[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==2))]<-"ACC"
#----------------------------------------------------
# A in C1 or G2
#----------------------------------------------------

mut.cols.1<-c("C1","G1","C2","G2")
mut.cols.2<-c("R1","Y3","R4","Y6")

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==2)
             &(rowSums(class.matrix[,c(mut.cols.1[c(1,4)])])==1)
)]<-1

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==2)
             &(rowSums(class.matrix[,c(mut.cols.1[c(1,4)])])==1)
             &(class.matrix[,"spacer"]==0)
)]<-3

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==2)
             &(rowSums(class.matrix[,c(mut.cols.1[c(1,4)])])==1)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
)]<-4

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==2)
             &(rowSums(class.matrix[,c(mut.cols.1[c(1,4)])])==1)
             &(class.matrix[,"spacer"]%in%c(1,2))
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
)]<-2

#----------------------------------------------------
# A in G1 or C2
#----------------------------------------------------

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==2)
             &(rowSums(class.matrix[,c(mut.cols.1[c(1,4)])])==0)
)]<-1

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==2)
             &(rowSums(class.matrix[,c(mut.cols.1[c(1,4)])])==0)
             &(class.matrix[,"spacer"]==0)
)]<-2

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==2)
             &(rowSums(class.matrix[,c(mut.cols.1[c(1,4)])])==0)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
)]<-3

#------------------------------------------------------

# 17) BCC

#----------------------------------------------------
mut.cols.1<-c("W1","W2","W3","W4")
mut.cols.2<-c("R1","Y3","R4","Y6")
labels[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==2))]<-"BCC"

#----------------------------------------------------
# B in W1 or W4
#----------------------------------------------------

mut.cols.1<-c("W1","W2","W3","W4")
mut.cols.2<-c("R1","Y3","R4","Y6")

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==2)
             &(rowSums(class.matrix[,c(mut.cols.1[c(1,4)])])==1)
)]<-1

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==2)
             &(rowSums(class.matrix[,c(mut.cols.1[c(1,4)])])==1)
             &(class.matrix[,"spacer"]==0)
)]<-3

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==2)
             &(rowSums(class.matrix[,c(mut.cols.1[c(1,4)])])==1)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
)]<-4

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==2)
             &(rowSums(class.matrix[,c(mut.cols.1[c(1,4)])])==1)
             &(class.matrix[,"spacer"]%in%c(1,2))
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
)]<-2


#----------------------------------------------------
# B in W2 or W3
#----------------------------------------------------


grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==2)
             &(rowSums(class.matrix[,c(mut.cols.1[c(1,4)])])==0)
)]<-1

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==2)
             &(rowSums(class.matrix[,c(mut.cols.1[c(1,4)])])==0)
             &(class.matrix[,"spacer"]==0)
)]<-2

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==2)
             &(rowSums(class.matrix[,c(mut.cols.1[c(1,4)])])==0)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
)]<-3



#--------------------------------------------------
### 18) CCC
#----------------------------------------------------

mut.cols<-c("R1","Y3","R4","Y6")
labels[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols])==3))]<-"CCC"
#----------------------------------------------------

mut.cols<-c("R1","Y3","R4","Y6")

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols])==3))]<-1

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols])==3)
             &(class.matrix[,"spacer"]==0))]<-3

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols])==3)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]==0))]<-4

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols])==3)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]%in%c(1,2)))]<-2

#--------------------------------------------------




#------------------------------------------------------

# 19) CCD

#----------------------------------------------------
mut.cols.1<-c("R2","R3","Y1","Y2","R5","R6","Y4","Y5")
mut.cols.2<-c("R1","Y3","R4","Y6")
labels[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==2))]<-"CCD"

#----------------------------------------------------
# D in first or last quarter
#----------------------------------------------------

mut.cols.1<-c("R2","R3","Y1","Y2","R5","R6","Y4","Y5")
mut.cols.2<-c("R1","Y3","R4","Y6")

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==2)
             &(rowSums(class.matrix[,c(mut.cols.1[c(1,2,7,8)])])==1)
)]<-1

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==2)
             &(rowSums(class.matrix[,c(mut.cols.1[c(1,2,7,8)])])==1)
             &(class.matrix[,"spacer"]==0)
)]<-3

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==2)
             &(rowSums(class.matrix[,c(mut.cols.1[c(1,2,7,8)])])==1)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
)]<-4

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==2)
             &(rowSums(class.matrix[,c(mut.cols.1[c(1,2,7,8)])])==1)
             &(class.matrix[,"spacer"]%in%c(1,2))
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
)]<-2


#----------------------------------------------------
# D in W2 or W3
#----------------------------------------------------


grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==2)
             &(rowSums(class.matrix[,c(mut.cols.1[c(1,2,7,8)])])==0)
)]<-1

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==2)
             &(rowSums(class.matrix[,c(mut.cols.1[c(1,2,7,8)])])==0)
             &(class.matrix[,"spacer"]==0)
)]<-2

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==1)
             &(rowSums(class.matrix[,mut.cols.2])==2)
             &(rowSums(class.matrix[,c(mut.cols.1[c(1,2,7,8)])])==0)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
)]<-3



#--------------------------------------------------
### 20) CCCC
#----------------------------------------------------

mut.cols<-c("R1","Y3","R4","Y6")
labels[which((class.matrix[,"n.mut.tot"]==4)
             &(rowSums(class.matrix[,mut.cols])==4))]<-"CCCC"
#----------------------------------------------------

mut.cols<-c("R1","Y3","R4","Y6")

grades[which((class.matrix[,"n.mut.tot"]==4)
             &(rowSums(class.matrix[,mut.cols])==4))]<-1

grades[which((class.matrix[,"n.mut.tot"]==4)
             &(rowSums(class.matrix[,mut.cols])==4)
             &(class.matrix[,"spacer"]==0))]<-2

grades[which((class.matrix[,"n.mut.tot"]==4)
             &(rowSums(class.matrix[,mut.cols])==4)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]==0))]<-3

grades[which((class.matrix[,"n.mut.tot"]==4)
             &(rowSums(class.matrix[,mut.cols])==4)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
             &(class.matrix[,"spacer"]%in%c(1,2)))]<-2

#-------------------------------------------------

#------------------------------------------------
### 21) ABC 

#First quarter
#---------------------------------------------------------

mut.cols.1<-c("C1","W1")
mut.cols.2<-c("R1","Y6")
labels[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1))]<-"ABC"
#-----------------------------------------------------------

grades[which((class.matrix[,"n.mut.tot"]==3)
       &(rowSums(class.matrix[,mut.cols.1])==2)
       &(rowSums(class.matrix[,mut.cols.2])==1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-3

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(class.matrix[,"spacer"]%in%c(1,2))
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-2

#-----------------------------------------------------

#Fourth quarter
#---------------------------------------------------------

mut.cols.1<-c("G2","W4")
mut.cols.2<-c("R1","Y6")
labels[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1))]<-"ABC"
#-----------------------------------------------------------

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-3

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(class.matrix[,"spacer"]%in%c(1,2))
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-2

#-----------------------------------------------------

#------------------------------------------------
### 22) BCD 

#First quarter
#---------------------------------------------------------

mut.cols.1<-c("W1","W1")
mut.cols.2<-c("R2","R3")
mut.cols.3<-c("R1","Y6")
labels[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,mut.cols.3])==1))]<-"BCD"
#-----------------------------------------------------------

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,mut.cols.3])==1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,mut.cols.3])==1)
             &(class.matrix[,"spacer"]==0))]<-2

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,mut.cols.3])==1)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-3

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,mut.cols.3])==1)
             &(class.matrix[,"spacer"]%in%c(1,2))
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-2

#--------------------------------------------------------
#Fourth quarter
#---------------------------------------------------------

mut.cols.1<-c("W4","W4")
mut.cols.2<-c("Y4","Y5")
mut.cols.3<-c("R1","Y6")
labels[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,mut.cols.3])==1))]<-"BCD"
#-----------------------------------------------------------

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,mut.cols.3])==1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,mut.cols.3])==1)
             &(class.matrix[,"spacer"]==0))]<-2

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,mut.cols.3])==1)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-3

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,mut.cols.3])==1)
             &(class.matrix[,"spacer"]%in%c(1,2))
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-2

#------------------------------------------------
#------------------------------------------------
### 23) ACD 

#First quarter
#--------------------------------------------------------------

mut.cols.1<-c("C1","C1")
mut.cols.2<-c("R2","R3")
mut.cols.3<-c("R1","Y6")
labels[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,mut.cols.3])==1))]<-"ACD"
#--------------------------------------------------------------

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,mut.cols.3])==1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,mut.cols.3])==1)
             &(class.matrix[,"spacer"]==0))]<-2

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,mut.cols.3])==1)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-3

#--------------------------------------------------------

#Fourth quarter
#--------------------------------------------------------------

mut.cols.1<-c("G2","G2")
mut.cols.2<-c("Y4","Y5")
mut.cols.3<-c("R1","Y6")
labels[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,mut.cols.3])==1))]<-"ACD"
#--------------------------------------------------------------

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,mut.cols.3])==1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,mut.cols.3])==1)
             &(class.matrix[,"spacer"]==0))]<-2

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(rowSums(class.matrix[,mut.cols.3])==1)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-3

#--------------------------------------------------------

#-------------------
### 24) CDD 
#-------------------

# all D and C in the 1st or 4th quarter

mut.cols.1<-c("R2","R3","Y4","Y5")
mut.cols.2<-c("R1","Y6")

labels[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1))]<-"CDD"


grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-2

#First quarter
#---------------------------------------------------------

mut.cols.1<-c("R2","R3")
mut.cols.2<-c("R1","Y6")
labels[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1))]<-"CDD"
#-----------------------------------------------------------

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(class.matrix[,"spacer"]==0))]<-2

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-3

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(class.matrix[,"spacer"]%in%c(1,2))
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-2

#-----------------------------------------------------

#Fourth quarter
#---------------------------------------------------------

mut.cols.1<-c("Y4","Y5")
mut.cols.2<-c("R1","Y6")
labels[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1))]<-"CDD"
#-----------------------------------------------------------

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1))]<-1

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(class.matrix[,"spacer"]==0))]<-2

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(class.matrix[,"spacer"]==0)
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-3

grades[which((class.matrix[,"n.mut.tot"]==3)
             &(rowSums(class.matrix[,mut.cols.1])==2)
             &(rowSums(class.matrix[,mut.cols.2])==1)
             &(class.matrix[,"spacer"]%in%c(1,2))
             &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-2

#-----------------------------------------------------

resume.matrix<-cbind(class.matrix,labels,grades)
if(length(which(resume.matrix[,"grades"]>0))>0){
resume.matrix.full<-resume.matrix[which(resume.matrix[,"grades"]>0),c("start","stop","spacer","n.mut.tot","sequence","mutations","WW1","WW2","labels","grades")]
}
}

pre.complete<-rbind(resume.matrix.full,resume.matrix.halves)

return(pre.complete)

}

