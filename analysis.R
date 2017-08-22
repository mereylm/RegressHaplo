library(readr)
library(xlsx)

P <- read.csv("~/Dropbox/Georgetown - Math & Stats/Haplotype Reconstruction/P.csv",header=FALSE)
P <- as.matrix(P)
y <- read_csv("~/Dropbox/Georgetown - Math & Stats/Haplotype Reconstruction/y.csv",col_names=FALSE)
y <- as.matrix(y)

rho_list <- c(rho_001=0.001,rho_01=0.01,rho_1=0.1,rho1=1)
full_output <- lapply(rho_list,analyze,y=y,P=P)

analyze <- function(rho,y,P) {
  inner_exit_status <- function(x) {
    #input is a row representing all inner loop zero/non-zero status
    #of a single element of the pi vector
    x <- x[!is.na(x)]
    hit_zero <- any(x==0)
    end_nonzero <- x[length(x)] == 1
    return(c(hit_zero,end_nonzero))
  }
  
  summarize_exit <- function(x) {
    #input is a matrix representing with 2 boolean columns:
    #Col (1): did that element ever hit zero in an inner loop of a single outer loop?
    #Col (2): did that element exit the single outer loop with non-zero status?
    exit_nonzero <- sum(x[,2])
    exit_nonzero_after_zero <- sum(x[,1] & x[,2])
    proportion <- exit_nonzero_after_zero/exit_nonzero
    return(c(exit_nonzero=exit_nonzero,proportion=proportion))
  }
  
  output <- penalized_regression.RegressHaplo(y,P,rho)
  # pi_record has dim=c(pi_len,inner_loop,stop_outer)
  # meaning that each matrix layer (i.e. [,,i]) represents
  # one outer loop iteration. In that layer, the column j is
  # the frequency (pi) vector returned during the jth iteration of
  # the inner loop
  
  # We want to know:
  # (1) when an element of pi is 0 during one inner loop iteration,
  #     how often does it come back to a positive value
  # (2) when an element of pi is 0 during the _last_ inner loop iteration of
  #     outer loop iteration i, how often will it be positive at the _last_
  #     inner loop iteration of outer loop iteration i+1
  
  pi_record <- output$pi_record
  pi_record_bin <- ifelse(pi_record>0,1,0)
  
  # (1)
  # find when an element of pi drops to 0 or comes back
  # if 1: element came back from being 0
  # if -1: element dropped to 0
  # if 0: element stayed either nonzero or zero
  inner_diff <- apply(pi_record_bin,c(1,3),diff)
  inner_diff <- aperm(inner_diff,c(2,1,3)) #rotate the forced output of 'apply'
  
  # Separate out the instances of dropping to zero vs coming back from zero
  drops <- apply(inner_diff,c(1,3), function(x) sum(x[x<0],na.rm=TRUE))
  come_backs <- apply(inner_diff,c(1,3), function(x) sum(x[x>0],na.rm=TRUE))
  
  # how often do elements that hit zero end up being non-zero?
  # let's start with the inner loop...that is, does it exit the inner loop as non-zero
  hit_zero_inner <- apply(pi_record_bin,c(1,3),inner_exit_status)
  hit_zero_inner <- aperm(hit_zero_inner,c(2,1,3)) #rotate the forced output of 'apply'
  
  # how many outer loops had the element hit zero and exit non-zero
  inner_exit_summary <- apply(hit_zero_inner, 3, summarize_exit)
  
  return(list(reconst=output,bool_mat=pi_record_bin,exit_summary=inner_exit_summary))
}




# Compute the total number of "flips", either dropping to zero or returning
changes <- apply(inner_diff,c(1,3),function(x) sum(abs(x),na.rm=TRUE))
flipping <- (changes-1)*(changes>0)
each_elem <- apply(flipping,1,sum)
all_flips <- sum(flipping)
sum(pi_record_bin,na.rm=TRUE)





# try different rho values, print to xlsx
write.xlsx(t(inner_exit_summary), file = "summary_results.xlsx",
           sheetName = "rho_1",append=TRUE)

# we would like to do cross-validation on the rho -- throw out rows of the P and y, train on the other ones and test on the hold out

