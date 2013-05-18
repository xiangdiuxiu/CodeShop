
myEval <- function(y, n, alpha){
  return(sum(as.matrix(dist(y))^alpha)/(2*n))
}


myImprove <- function(y, x, ny, nvy, alpha){
  tempdata <- cbind(y,x)
  tempdata <- tempdata[order(tempdata[,ncol(tempdata)]),] 
  ytemp <- tempdata[,1:nvy]
  xtemp <- tempdata[,nvy+1]
  dy <- as.matrix(dist(y))^alpha
  best <- 0
  which <- 0
  for(i in 1:(ny-1)){
    left_sum <- sum(dy[(1:i),(1:i)])/(2*i)
    right_sum <- sum(dy[-(1:i),-(1:i)])/(2*(ny-i))
    SUM <- left_sum+right_sum
    if(is.na(SUM)) {show(i)}
    if(is.nan(SUM)) {show(i)}
    if(SUM>best && i>5 && (ny-i)>5){
      which <- i
      best <- SUM
    }
  }
  
  if(best>0){
    return(list(bv=best,sv=xtemp[which]))
  }
  
}

chooseSplit <- function(y, x, ny, nvy, nvx, alpha){
  best <- 0
  bestsp <- 0
  bestspv <- 0
  for (i in 1:nvx) {
    res <- myImprove(y, x[i], ny, nvy, alpha)
    if(res$best>best){
      best <- res$best
      bestsp <- i
      bestspv <- res$sv
    }
  }
  return(list(bv=best, sp=bestsp, spv=bestspv))
}

myRpart <- function(data, alpha, nvy, nvx, nodeNumber, FnodeNumber, min, sp, spvalue){
  ytemp <- data[, 1:nvy]
  xtemp <- data[, (nvy+1):(nvy+nvx)]
  ny <- nrow(ytemp) # number of sample
  if(nodeNumber==0){
    nodeType = "Root"
  }
  
  if(ny<min){
    nodeType = "Leaf"
    show(list(NODE=nodeNumber, FNODE=FnodeNumber, SplitVal=sp, Splitvalue=spvalue,Nsample=ny, value=myEval(ytemp, ny, alpha), TYPE=nodeType))
    return()
  }
  #print current node to screen
  show(list(NODE=nodeNumber, FNODE=FnodeNumber, SplitVal=sp, Splitvalue=spvalue,Nsample=ny, value=myEval(ytemp, ny, alpha), TYPE=nodeType))
  split <- chooseSplit(ytemp, xtemp, ny, nvy, nvx, alpha)
  spval <- split$sp
  spvalue <- split$spv
  myRpart(data[data[,nvy+spval]<=spvalue,], alpha, nvy, nvx, 2*(nodeNumber)+1, nodeNumber, min, spval, spvalue)
  myRpart(data[data[,nvy+spval]>spvalue,], alpha, nvy, nvx, 2*(nodeNumber)+2, nodeNumber, min, spval, spvalue)
  
}
