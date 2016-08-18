#' @title Function p53.track
#' @param seq A character string containing the sequence. The sequence must be composed exclusively of DNA bases (a,c,t,g)
#' @param seqname A character string containing the name of the sequence. The default is an empty string.
#' @param plot A boolean value indicating whether the results should be displayed on a plot. The default is FALSE.
#' @return A dataframe containing one row for each responsive elements located on the input sequence. For each element, the following infomation is provided: 
#' #' \itemize{
#' \item ID: sequence ID (e.g. CDKN1A)
#' \item start: coordinate of the RE's start with respect to the input sequence (e.g. 246)
#' \item stop: relative coordinate of the RE's stop with respect to the input sequence (e.g. 266)
#' \item spacer: spacer length of the RE (e.g. 10)
#' \item n.mismatches: number of mismatches with respect to the canonical p53 RE (e.g. 2)
#' \item sequence: RE sequence (including the spacer, e.g. gaacctgcttcatttaaatggagcatgtgt)
#' \item mismatch.string: string encoding the position of mismatches in the RE with respect to the canonical one (0 = match, 1 = mismatch, n= spacer, e.g. 0000100000nnnnnnnnnn0000000010)
#' \item WW1: sequence of WW1 in the RE (e.g. CT)
#' \item WW2: sequence of WW2 in the RE (e.g. AT)
#' \item label: mismatch label assigned to the RE by p53retriever, according to the number, the nature and the position of mismatches (0 = no mismatches, A = mismatch in C or G positions, B = mismatch in W positions, C = mismatch in R or Y "external" positions, D = mismatch in R or Y "internal" positions, 3Q = three-quarter site, half = half site)
#' \item grade: functional score assigned to the RE by p53retriever, according to the number, the nature and the position of mismatches (5 = high, 4 = moderate, 3 = slight, 2 = poor, 1= unlikely)
#' }
#' @export
#' @import Biostrings
#' @description This function locates candidate p53 responsive elements on a DNA sequence. The rules formalized in this model are manually curated and based on observations obtained during several years of experiments using the yeast-based assays, including the results presented in (Inga et al, MCB 2002; Tomso et al, PNAS 2005; Jegga et al, PNAS 2008; Jordan et al, PloS Genetics 2008; Menendez et al, PNAS 2010).
#' @examples
#' data(CDKN1A)
#' p53track(CDKN1A,seqname="CDKN1A",plot=TRUE)

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
