#' Return the maximum scores using jaspar p53 PWMs (pwm.1 derived from MA0106.1, built on SELEX data, pwm.2 derived from MA0106.2, built on ChIPseq data). The score ranges from 0 to 1
#' 
#' @param seq A character string containing the sequence. The sequence must be composed exclusively of DNA bases (a,c,t,g)
#' @return A vector containing the two maximum scores given by the two jaspar PWMs.


#--------------------------------------------------------------------------
# Return the maximum score using jaspar PWMs. The score ranges from 0 to 1
#--------------------------------------------------------------------------

p53maxJaspar<-function(seq){
  
  data(p53.jaspar.pwm)
  subject<-DNAString(seq)
  
  j1.hits <- matchPWM(j1.pwm, subject, min.score=0,with.score=TRUE)
  if(length(j1.hits)>0) {
    max.j1.score<-max(mcols(j1.hits)$score)
  } else {
    max.j1.score<-0
  }  
  
  j2.hits <- matchPWM(j2.pwm, subject, min.score=0,with.score=TRUE)
  if(length(j2.hits)>0) {
    max.j2.score<-max(mcols(j2.hits)$score)
  } else {
    max.j2.score<-0
  }
  scores<-c(max.j1.score,max.j2.score)
  names(scores)<-c("jaspar.selex","jaspar.chipseq")
  return (scores)
}
