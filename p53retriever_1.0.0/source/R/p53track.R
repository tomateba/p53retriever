#' @title Function p53.track
#' @param seq A character string containing the sequence. The sequence must be composed exclusively of DNA bases (a,c,t,g)
#' @param seqname A character string containing the name of the sequence. The default is an empty string.
#' @param plot A boolean value indicating whether the results should be displayed on a plot. The default is FALSE.
#' @return A dataframe containing the responsive elements located on the input sequence. For each element, this infomation is provided: start and stop coordinates, spacer length, number of mismatches, sequence, mismatch string, mismatch label and grade. If plot is set to TRUE, a plot displaying the hits is also created.
#' @export
#' @import Biostrings
#' @description This function locates candidate p53 responsive elements on a DNA sequence. The rules formalized in this model are manually curated and based on observations obtained during several years of experiments using the yeast-based assays, including the results presented in (Inga et al, MCB 2002; Tomso et al, PNAS 2005; Jegga et al, PNAS 2008; Jordan et al, PloS Genetics 2008; Menendez et al, PNAS 2010).
#' @examples
#' data(CDKN1A)
#' p53track(CDKN1A.promoter.seq,seqname="CDKN1A",plot=TRUE)

p53track<-function(seq,seqname="",plot=F){

seq.ini<-seq
if((all(unique(unlist(strsplit(as.character(seq.ini),split="")))%in%c("a","c","g","t","A","C","G","T")))==FALSE){
  stop('The sequence contains invalid characters (only DNA bases are allowed)')
}

full<-p53sf(seq.ini) 
q3<-p53s3q(seq.ini)

complete<-rbind(full,q3)
  
#-----------------------------------------
# Clean data
#-----------------------------------------
  
p53.clean<-NULL
  
if(class(complete)!="NULL"){
  
  pre.ord.matrix<-complete[order(complete[,"start"]),]
  ord.matrix<-pre.ord.matrix
  
  if(nrow(pre.ord.matrix)>1){
    i<-2
    while(i<=nrow(ord.matrix)){
      if(ord.matrix[i,"start"]<ord.matrix[i-1,"stop"]){
        if(ord.matrix[i,"grades"]>ord.matrix[i-1,"grades"]) {ord.matrix<-ord.matrix[-(i-1),]}
        else {if(ord.matrix[i,"grades"]<ord.matrix[i-1,"grades"]) {ord.matrix<-ord.matrix[-(i),]}
              else {
                if(nchar(as.character(ord.matrix[i,"sequence"]))>nchar(as.character(ord.matrix[i-1,"sequence"]))){ord.matrix<-ord.matrix[-(i-1),]}
                else{ord.matrix<-ord.matrix[-(i),]}}}}
      else {i<-i+1}
    }
  }
  rownames(ord.matrix)<-NULL
  ID<-as.character(rep(seqname,nrow(ord.matrix)))
  p53.clean<-cbind(ID,ord.matrix)
  colnames(p53.clean)<-c("ID","start","stop","spacer","n.mismatches","sequence","mismatch.string","WW1","WW2","label","grade")
} 
  
if (plot==T){
  p53plot(p53.clean,0,nchar(seq.ini))  
}

  return(p53.clean) 
}
