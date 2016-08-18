#' Locate candidate p53 responsive elements (3/4 sites) on a DNA sequence
#' 
#' @param seq.ini A character string containing the sequence. The sequence must be composed exclusively of DNA bases (a,c,t,g)
#' @return A dataframe containing the responsive elements located on the input sequence.
#' @export 

p53s3q<-function(seq.ini){

data(bkg.models)
subject<-DNAString(seq.ini)

#-------------------------------------------------
### 3/4 123Q
#-------------------------------------------------

matrix.123<-resume.matrix.123<-NULL
matches.123<- matchPWM(pwm.2, subject, min.score="92%")

if(length(matches.123)>0){
  
  sequence<-as.character(matches.123)
  start.seq<-start(matches.123)
  stop.seq<-end(matches.123)
  ww<-unlist(lapply(sequence, function(x) substr(x,5,6)))
  matrix.123<-cbind(sequence,start.seq,stop.seq,ww)
  colnames(matrix.123)<-c("seq","start","stop","WW")
  
  work.matrix<-matrix.123
  pattern.matrix<-matrix(0,nrow=nrow(work.matrix),ncol=15)
  ww.matrix<-matrix("",nrow=nrow(work.matrix),ncol=2)
  colnames(pattern.matrix)<-c("R1","R2","R3","C1","W1","W2","G1","Y1","Y2","Y3",
                           "R4","R5","R6","C2","W3")
  colnames(ww.matrix)<-c("WW1","WW2")
  n.mut.tot<-vector(mode="numeric",length=nrow(work.matrix))
  spacer<-rep(0,nrow(work.matrix))
  sequence<-work.matrix[,"seq"]
  double.pwm<-pwm.2
  
  for(p in (1:nrow(work.matrix))){
    pair.seq<-work.matrix[p,"seq"]
    line<-vector(mode="numeric",length=15)
    for (l in (1:15)){
      letter<-substr(pair.seq,l,l)
      line[l]<-as.numeric(double.pwm[as.character(letter),l]=="0")
    }
    pattern.matrix[p,]<-line
    ww.matrix[p,1]<-substr(pair.seq,5,6)
    ww.matrix[p,2]<-""
    n.mut.tot[p]<-sum(line[1:15]==1)
  }
  
  pattern<-as.data.frame(pattern.matrix)
  ww<-as.data.frame(ww.matrix)
  
  
  mutations<-apply(pattern,1,function(x) paste(paste(x[1:15],collapse=""),"nnnnn",sep=""))
  
  smart.matrix<-cbind(matrix(work.matrix[,c("start","stop")],nrow=nrow(work.matrix)),spacer,n.mut.tot,sequence,mutations,ww,pattern)
  colnames(smart.matrix)[1:2]<-c("start","stop")
  
  smart.matrix[,"start"]<-as.numeric(work.matrix[,"start"])
  smart.matrix[,"stop"]<-as.numeric(work.matrix[,"stop"])
  
  class.matrix<-transform(smart.matrix,
                          sequence=as.character(sequence),
                          mutations=as.character(mutations),
                          WW1=as.character(WW1),
                          WW2=as.character(WW2))
  
  grades<-rep(0,length=nrow(class.matrix))
  labels<-rep("",length=nrow(class.matrix))
  
  
  #  0 mutations
  #-----------------------------------------------------
  labels[which(class.matrix[,"n.mut.tot"]==0)]<-"3Q.I-II-III"
  #-----------------------------------------------------
  
  grades[which(class.matrix[,"n.mut.tot"]==0)]<-3
  grades[which((class.matrix[,"n.mut.tot"]==0)
               &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-4
  
  #1) C mutation
  #----------------------------------------------------
  
  labels[which((class.matrix[,"n.mut.tot"]==1)
               &(class.matrix[,"R1"]==1|class.matrix[,"Y3"]==1|class.matrix[,"R4"]==1)
  )]<-"3Q.I-II-III.C"
  #----------------------------------------------------
  
  #1.1) C in the first quarter
  
  grades[which((class.matrix[,"n.mut.tot"]==1)
               &(class.matrix[,"R1"]==1|class.matrix[,"Y3"]==1|class.matrix[,"R4"]==1)
               &(class.matrix[,"R1"]==1)
  )]<-3
  
  grades[which((class.matrix[,"n.mut.tot"]==1)
               &(class.matrix[,"R1"]==1|class.matrix[,"Y3"]==1|class.matrix[,"R4"]==1)
               &(class.matrix[,"R1"]==1)
               &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
  )]<-4
  
  #1.2) C in the second or third quarter
  
  grades[which((class.matrix[,"n.mut.tot"]==1)
               &(class.matrix[,"R1"]==1|class.matrix[,"Y3"]==1|class.matrix[,"R4"]==1)
               &(class.matrix[,"R1"]==0)
  )]<-2
  
  grades[which((class.matrix[,"n.mut.tot"]==1)
               &(class.matrix[,"R1"]==1|class.matrix[,"Y3"]==1|class.matrix[,"R4"]==1)
               &(class.matrix[,"R1"]==0)
               &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
  )]<-3
  
  
  resume.matrix<-cbind(class.matrix,labels,grades)
  
  if(length(which(resume.matrix[,"grades"]>0))>0){
    resume.matrix.123<-resume.matrix[which(resume.matrix[,"grades"]>0),c("start","stop","spacer","n.mut.tot","sequence","mutations","WW1","WW2","labels","grades")]
  }
}

#-------------------------------------------------
### 3/4 234Q
#-------------------------------------------------

matrix.234<-resume.matrix.234<-NULL
matches.234<- matchPWM(pwm.3, subject, min.score="92%")

if(length(matches.234)>0){
  
  sequence<-as.character(matches.234)
  start.seq<-start(matches.234)
  stop.seq<-end(matches.234)
  ww<-unlist(lapply(sequence, function(x) substr(x,10,11)))
  matrix.234<-cbind(sequence,start.seq,stop.seq,ww)
  colnames(matrix.234)<-c("seq","start","stop","WW")
  
  work.matrix<-matrix.234
  pattern.matrix<-matrix(0,nrow=nrow(work.matrix),ncol=15)
  ww.matrix<-matrix("",nrow=nrow(work.matrix),ncol=2)
  colnames(pattern.matrix)<-c("W2","G1","Y1","Y2","Y3","R4","R5","R6","C2","W3","W4","G2","Y4","Y5","Y6")
  colnames(ww.matrix)<-c("WW1","WW2")
  n.mut.tot<-vector(mode="numeric",length=nrow(work.matrix))
  spacer<-rep(0,nrow(work.matrix))
  sequence<-work.matrix[,"seq"]
  double.pwm<-pwm.3
  
  for(p in (1:nrow(work.matrix))){
    pair.seq<-work.matrix[p,"seq"]
    line<-vector(mode="numeric",length=15)
    for (l in (1:15)){
      letter<-substr(pair.seq,l,l)
      line[l]<-as.numeric(double.pwm[as.character(letter),l]=="0")
    }
    pattern.matrix[p,]<-line
    ww.matrix[p,1]<-""
    ww.matrix[p,2]<-substr(pair.seq,10,11)
    n.mut.tot[p]<-sum(line[1:15]==1)
  }
  
  pattern<-as.data.frame(pattern.matrix)
  ww<-as.data.frame(ww.matrix)
  
  
  mutations<-apply(pattern,1,function(x) paste("nnnnn",paste(x[1:15],collapse=""),sep=""))
  
  smart.matrix<-cbind(matrix(work.matrix[,c("start","stop")],nrow=nrow(work.matrix)),spacer,n.mut.tot,sequence,mutations,ww,pattern)
  colnames(smart.matrix)[1:2]<-c("start","stop")
  
  smart.matrix[,"start"]<-as.numeric(work.matrix[,"start"])
  smart.matrix[,"stop"]<-as.numeric(work.matrix[,"stop"])
  
  class.matrix<-transform(smart.matrix,
                          sequence=as.character(sequence),
                          mutations=as.character(mutations),
                          WW1=as.character(WW1),
                          WW2=as.character(WW2))
  
  grades<-rep(0,length=nrow(class.matrix))
  labels<-rep("",length=nrow(class.matrix))
  
  #  0 mutations 
  #-----------------------------------------------------
  labels[which(class.matrix[,"n.mut.tot"]==0)]<-"3Q.II-III-IV"
  #-----------------------------------------------------
  
  grades[which(class.matrix[,"n.mut.tot"]==0)]<-3
  
  grades[which((class.matrix[,"n.mut.tot"]==0)
               &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-4
  
  # 1 C mutation
  #----------------------------------------------------
  
  labels[which((class.matrix[,"n.mut.tot"]==1)
               &(class.matrix[,"Y3"]==1|class.matrix[,"R4"]==1|class.matrix[,"Y6"]==1)
  )]<-"3Q.II-III-IV.C"
  #----------------------------------------------------
  
  # C in the last quarter
  
  grades[which((class.matrix[,"n.mut.tot"]==1)
               &(class.matrix[,"Y3"]==1|class.matrix[,"R4"]==1|class.matrix[,"Y6"]==1)
               &(class.matrix[,"Y6"]==1)
  )]<-3
  
  grades[which((class.matrix[,"n.mut.tot"]==1)
               &(class.matrix[,"Y3"]==1|class.matrix[,"R4"]==1|class.matrix[,"Y6"]==1)
               &(class.matrix[,"Y6"]==1)
               &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
  )]<-4
  
  # C in the second or third quarter
  
  grades[which((class.matrix[,"n.mut.tot"]==1)
               &(class.matrix[,"Y3"]==1|class.matrix[,"R4"]==1|class.matrix[,"Y6"]==1)
               &(class.matrix[,"Y6"]==0)
  )]<-2
  
  grades[which((class.matrix[,"n.mut.tot"]==1)
               &(class.matrix[,"Y3"]==1|class.matrix[,"R4"]==1|class.matrix[,"Y6"]==1)
               &(class.matrix[,"Y6"]==0)
               &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
  )]<-3
  
  resume.matrix<-cbind(class.matrix,labels,grades)
  
  if(length(which(resume.matrix[,"grades"]>0))>0){
    resume.matrix.234<-resume.matrix[which(resume.matrix[,"grades"]>0),c("start","stop","spacer","n.mut.tot","sequence","mutations","WW1","WW2","labels","grades")]
  }
  
}


#-------------------------------------------------
### 3/4 134Q
#-------------------------------------------------

matrix.134<-resume.matrix.134<-NULL
matches.134<- matchPWM(pwm.4, subject, min.score="92%")


if(length(matches.134)>0){
  
  sequence<-as.character(matches.134)
  start.seq<-start(matches.134)
  stop.seq<-end(matches.134)
  ww<-unlist(lapply(sequence, function(x) substr(x,15,16)))
  matrix.134<-cbind(sequence,start.seq,stop.seq,ww)
  colnames(matrix.134)<-c("seq","start","stop","WW")
  
  work.matrix<-matrix.134
  pattern.matrix<-matrix(0,nrow=nrow(work.matrix),ncol=20)
  ww.matrix<-matrix("",nrow=nrow(work.matrix),ncol=2)
  colnames(pattern.matrix)<-c("R1","R2","R3","C1","W1","W2","G1","Y1","Y2","Y3",
                              "R4","R5","R6","C2","W3","W4","G2","Y4","Y5","Y6")
  colnames(ww.matrix)<-c("WW1","WW2")
  n.mut.tot<-vector(mode="numeric",length=nrow(work.matrix))
  spacer<-rep(0,nrow(work.matrix))
  sequence<-work.matrix[,"seq"]
  double.pwm<-pwm.4
  
  for(p in (1:nrow(work.matrix))){
    pair.seq<-work.matrix[p,"seq"]
    line<-vector(mode="numeric",length=20)
    for (l in (1:20)){
      letter<-substr(pair.seq,l,l)
      line[l]<-as.numeric(double.pwm[as.character(letter),l]=="0")
    }
    pattern.matrix[p,]<-line
    ww.matrix[p,1]<-""
    ww.matrix[p,2]<-substr(pair.seq,15,16)
    n.mut.tot[p]<-sum(line[c(1:5,11:20)]==1)
  }
  
  pattern<-as.data.frame(pattern.matrix)
  ww<-as.data.frame(ww.matrix)
  
  
  mutations<-apply(pattern,1,function(x) paste(paste(x[1:5],collapse=""),"nnnnn",paste(x[11:20],collapse=""),sep=""))
  
  smart.matrix<-cbind(matrix(work.matrix[,c("start","stop")],nrow=nrow(work.matrix)),spacer,n.mut.tot,sequence,mutations,ww,pattern)
  colnames(smart.matrix)[1:2]<-c("start","stop")
  
  smart.matrix[,"start"]<-as.numeric(work.matrix[,"start"])
  smart.matrix[,"stop"]<-as.numeric(work.matrix[,"stop"])
  
  class.matrix<-transform(smart.matrix,
                          sequence=as.character(sequence),
                          mutations=as.character(mutations),
                          WW1=as.character(WW1),
                          WW2=as.character(WW2))
  
  grades<-rep(0,length=nrow(class.matrix))
  labels<-rep("",length=nrow(class.matrix))
  
  #  0 mutations
  #-----------------------------------------------------
  labels[which(class.matrix[,"n.mut.tot"]==0)]<-"3Q.I-III-IV"
  #-----------------------------------------------------
  
  grades[which(class.matrix[,"n.mut.tot"]==0)]<-2
  
  grades[which((class.matrix[,"n.mut.tot"]==0)
               &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-3
  
  # 1 C mutation
  #----------------------------------------------------
  
  labels[which((class.matrix[,"n.mut.tot"]==1)
               &(class.matrix[,"R1"]==1|class.matrix[,"R4"]==1|class.matrix[,"Y6"]==1)
  )]<-"3Q.I-III-IV.C"
  #----------------------------------------------------
  
  grades[which((class.matrix[,"n.mut.tot"]==1)
               &(class.matrix[,"R1"]==1|class.matrix[,"R4"]==1|class.matrix[,"Y6"]==1)
  )]<-2
  
  grades[which((class.matrix[,"n.mut.tot"]==1)
               &(class.matrix[,"R1"]==1|class.matrix[,"R4"]==1|class.matrix[,"Y6"]==1)
               &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
  )]<-3
  
  
  resume.matrix<-cbind(class.matrix,labels,grades)
  
  if(length(which(resume.matrix[,"grades"]>0))>0){
    resume.matrix.134<-resume.matrix[which(resume.matrix[,"grades"]>0),c("start","stop","spacer","n.mut.tot","sequence","mutations","WW1","WW2","labels","grades")]
  }
  
}

#-------------------------------------------------
### 3/4 124Q
#-------------------------------------------------

matrix.124<-resume.matrix.124<-NULL
matches.124<- matchPWM(pwm.5, subject, min.score="92%")

if(length(matches.124)>0){
  
  sequence<-as.character(matches.124)
  start.seq<-start(matches.124)
  stop.seq<-end(matches.124)
  ww<-unlist(lapply(sequence, function(x) substr(x,5,6)))
  matrix.124<-cbind(sequence,start.seq,stop.seq,ww)
  colnames(matrix.124)<-c("seq","start","stop","WW")
  
  work.matrix<-matrix.124
  pattern.matrix<-matrix(0,nrow=nrow(work.matrix),ncol=20)
  ww.matrix<-matrix("",nrow=nrow(work.matrix),ncol=2)
  colnames(pattern.matrix)<-c("R1","R2","R3","C1","W1","W2","G1","Y1","Y2","Y3",
                              "R4","R5","R6","C2","W3","W4","G2","Y4","Y5","Y6")
  colnames(ww.matrix)<-c("WW1","WW2")
  n.mut.tot<-vector(mode="numeric",length=nrow(work.matrix))
  spacer<-rep(0,nrow(work.matrix))
  sequence<-work.matrix[,"seq"]
  double.pwm<-pwm.5
  
  for(p in (1:nrow(work.matrix))){
    pair.seq<-work.matrix[p,"seq"]
    line<-vector(mode="numeric",length=20)
    for (l in (1:20)){
      letter<-substr(pair.seq,l,l)
      line[l]<-as.numeric(double.pwm[as.character(letter),l]=="0")
    }
    pattern.matrix[p,]<-line
    ww.matrix[p,1]<-substr(pair.seq,5,6)
    ww.matrix[p,2]<-""
    n.mut.tot[p]<-sum(line[c(1:10,16:20)]==1)
  }
  
  pattern<-as.data.frame(pattern.matrix)
  ww<-as.data.frame(ww.matrix)
  
  mutations<-apply(pattern,1,function(x) paste(paste(x[1:10],collapse=""),"nnnnn",paste(x[16:20],collapse=""),sep=""))
  
  smart.matrix<-cbind(matrix(work.matrix[,c("start","stop")],nrow=nrow(work.matrix)),spacer,n.mut.tot,sequence,mutations,ww,pattern)
  colnames(smart.matrix)[1:2]<-c("start","stop")
  
  smart.matrix[,"start"]<-as.numeric(work.matrix[,"start"])
  smart.matrix[,"stop"]<-as.numeric(work.matrix[,"stop"])
  
  class.matrix<-transform(smart.matrix,
                          sequence=as.character(sequence),
                          mutations=as.character(mutations),
                          WW1=as.character(WW1),
                          WW2=as.character(WW2))
  
  grades<-rep(0,length=nrow(class.matrix))
  labels<-rep("",length=nrow(class.matrix))
  
  #  0 mutations
  #-----------------------------------------------------
  labels[which(class.matrix[,"n.mut.tot"]==0)]<-"3Q.I-II-IV"
  #-----------------------------------------------------
  
  grades[which(class.matrix[,"n.mut.tot"]==0)]<-2
  
  grades[which((class.matrix[,"n.mut.tot"]==0)
               &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT")))]<-3
  
  # 1 C mutation
  #----------------------------------------------------
  
  labels[which((class.matrix[,"n.mut.tot"]==1)
               &(class.matrix[,"R1"]==1|class.matrix[,"Y3"]==1|class.matrix[,"Y6"]==1)
  )]<-"3Q.I-II-IV.C"
  #----------------------------------------------------
  
  grades[which((class.matrix[,"n.mut.tot"]==1)
               &(class.matrix[,"R1"]==1|class.matrix[,"Y3"]==1|class.matrix[,"Y6"]==1)
  )]<-2
  
  grades[which((class.matrix[,"n.mut.tot"]==1)
               &(class.matrix[,"R1"]==1|class.matrix[,"Y3"]==1|class.matrix[,"Y6"]==1)
               &((class.matrix[,"WW1"]=="AT")|(class.matrix[,"WW2"]=="AT"))
  )]<-3
  
  
  resume.matrix<-cbind(class.matrix,labels,grades)
  
  if(length(which(resume.matrix[,"grades"]>0))>0){
    resume.matrix.124<-resume.matrix[which(resume.matrix[,"grades"]>0),c("start","stop","spacer","n.mut.tot","sequence","mutations","WW1","WW2","labels","grades")]
  }
  
}

#--------------------------------------------------

pre.complete<-rbind(resume.matrix.123,resume.matrix.234,resume.matrix.134,resume.matrix.124)

return(pre.complete)

}
