
round_it <- function( x )
{
  if(isTRUE((x - floor(x))>=0.5))
    ceiling(x)
  else
    floor(x)
}


nominal_to_binary <- function( data )
{
  result = NULL
  for (i in 1:ncol(data))
  {
     #print(i)
     if (is.numeric( data[,i] ) )
     {
        if (is.null(result))
          result = data.frame(data[,i])
        else
          result = data.frame(result, data[,i])
        colnames(result)[ncol(result)] <- colnames(data)[i]
     }
     else
     {
        vals = unique(data[,i])
        for (j in 1:length(vals))
        {
           #print(j)
           bins = c()
           for (k in 1:nrow(data))
           {
              if(data[,i][k] == vals[j])
                bins = c(bins,1)
              else
                bins = c(bins,0)
           }
           #print(bins)
           if (is.null(result))
             result = data.frame(bins)
           else
             result = data.frame(result, bins)
           colnames(result)[ncol(result)] <- paste(colnames(data)[i],"is",vals[j])
           if (length(vals)==2) break
        }
     }
  }
  #print(head(result))
  result
}

process_data <- function( data, colnames=NULL )
{
  data.num <- as.data.frame(data)
  if (!is.null(colnames))
  {
    data.num = subset(data.num, select = colnames)
  }
  if (!is.numeric(data.num))
  {
    data.num = nominal_to_binary(data.num)
  }
  if(any(is.na(data.num)))
  {
    require("gam")
   	data.repl = na.gam.replace(data.num)
  }
  else
  	data.repl = data.num
  data.repl
}

cluster <- function( data, min=10, max=15 )
{
  require("vegan")
  max <- min(max,nrow(unique(data)))
  max <- min(max,nrow(data)-1)
  if (min>max)
    min=max
  print(paste("cascade k-means ",min," - ",max))
  s = cascadeKM(data,min,max,iter=30)
  m = max.col(s$results)[2]
  print(paste("best k-means clustering result: ",((m-1)+min)," num clusters"))
  cbind(s$partition[,m])
}

stratified_split <- function( data, ratio=0.3, method="cluster", colnames=NULL )
{
    data.processed = as.matrix(process_data( data, colnames ))
    print(paste("split using #features: ",ncol(data.processed)))
    if (method == "samplecube")
    {
      require("sampling")
      # adjust ratio to make samplecube return exact number of samples
      ratio = round_it(nrow(data.processed)*ratio)/nrow(data.processed)
      pik = rep(ratio,times=nrow(data.processed))
      data.strat = cbind(pik,data.processed)
      samplecube(data.strat,pik,order=2,comment=F)
    }
    else if (method == "cluster")
    {
      cl = cluster(data.processed)
#      require("caret")
#      res = createDataPartition(cl,p=ratio)
#      split = rep(1, times=nrow(data))
#      for (j in 1:nrow(data))
#        if ( is.na(match(j,res$Resample1)) )
#          split[j]=0
#      split
      require("sampling")
      stratified_split(cl,ratio,"samplecube")
    }
    else
      stop("unknown method")
}

stratified_k_fold_split <- function( data, num_folds=10, method="cluster", colnames=NULL )
{
  print(paste(num_folds,"-fold-split, data-size",nrow(data)))
  data.processed = as.matrix(process_data( data, colnames ))
  print(paste("split using #features: ",ncol(data.processed)))
  if (method == "samplecube")
  {
    folds = rep(0, times=nrow(data))
    for (i in 1:(num_folds-1))
    {
      require("sampling")
      prop = 1/(num_folds-(i-1))
      print(paste("fold",i,"/",num_folds," prop",prop))
      pik = rep(prop,times=nrow(data))
      for (j in 1:nrow(data))
        if(folds[j]!=0)
          pik[j]=0
      data.strat = cbind(pik,data.processed)
      s<-samplecube(data.strat,pik,order=2,comment=F)
      print(paste("fold size: ",sum(s)))
      for (j in 1:nrow(data))
        if (s[j] == 1)
          folds[j]=i
    }
    for (j in 1:nrow(data))
      if (folds[j] == 0)
        folds[j]=num_folds
    folds
  }
  else if (method == "cluster")
  {
    require("TunePareto")
    cl = cluster(data.processed)
    res = generateCVRuns(cl,ntimes=1,nfold=num_folds)
    folds = rep(0, times=nrow(data))
    for (i in 1:num_folds)
      for(j in 1:length(res[[1]][[i]]))
        folds[res[[1]][[i]][j]]=i
    folds
  }
  else
    stop("unknown method")
}

duplicate_indices <- function( data ) {
  indices = 1:nrow(data) 
  z = data
  duplicate_index = anyDuplicated(z) 
  while(duplicate_index) {
    duplicate_to_index = anyDuplicated(z[1:duplicate_index,],fromLast=T)
    #print(paste(duplicate_index,'is dupl to',duplicate_to_index))
    indices[duplicate_index] <- duplicate_to_index
    z[duplicate_index,] <- paste('123$§%',duplicate_index)
    duplicate_index = anyDuplicated(z) 
  } 
  indices 
}

add_duplicates <- function( data, dup_indices ) { 
  result = data[1,]
  for(i in 2:length(dup_indices)) { 
    row = data[rownames(data)==dup_indices[i],]
    if(length(row)==0)
       stop(paste('index ',i,' dup-index ',dup_indices[i],'not found in data'))
    result = rbind(result, row) 
  } 
  rownames(result)<-NULL 
  result 
}

sammon_duplicates <- function( data, ... ) { 
  di <- duplicate_indices(data)
  print(di)
  u <- unique(data) 
  print(paste('unique data points',nrow(u),'of',nrow(data)))
  if(nrow(u) <= 4) stop("number of unqiue datapoints <= 4")
  points_unique <- sammon(dist(u), ...)$points
  if (nrow(u)<nrow(data))
  {
    points <- add_duplicates(points_unique, di) 
    points 
  }
  else
  {
    points_unique
  }
}

plot_pre_process <- function( data, method="pca" )
{
  data.processed = process_data( data )
  if (method == "pca")
  {
    data.pca <- prcomp(data.processed, scale=TRUE)
    as.data.frame(data.pca$x)[1:2]
  }
  else if (method == "smacof")
  {
    require("smacof")
    data.emb <- smacofSym(dist(data.processed, method = "euclidean"), ndim=2, verbose=T)
    data.emb$conf
  }
  else if (method == "sammon")
  {
    require("MASS")
    sammon_duplicates(data.processed, k=2)
  }
  else
    stop("unknown method")
}

plot_split <- function( data, split, names=NULL, ... )
{
  if (ncol(data)!=2 || !is.numeric(data[,1]) || !is.numeric(data[,2]))
    stop("data not suitable for plotting, plot_pre_process() first")

  plot( NULL, xlim = extendrange(data[,1]), ylim = extendrange(data[,2]), ... )
  if (is.null(names))
    names <- c("split 1","split 2")
  colos = as.double(rep(2:(max(split)+2)))
  legend("topleft",names,pch=2,col=colos)

  for (j in max(split):0)
  {
    set = c()
    for (i in 1:nrow(data))
      if (split[i] == j)
        set = c(set,i)
    points(data[set,], pch = 2, col=(j+2))
  }
}

#a<-matrix(rnorm(100, mean=50,  sd=4), ncol=5)
#b<-matrix(rnorm(5000, mean=0, sd=10), ncol=5)
#data<-rbind(a,b)
#c<-matrix(rnorm(50, mean=-50, sd=2), ncol=5)
#data<-rbind(data,c)
#data=iris
#split = stratified_k_fold_split(data, num_folds=3)
#split = stratified_split(data, ratio=0.33, method="cluster")
#print(sum(split))
#plot_split(plot_pre_process(data),split,c("training","test"))

#cl = cluster(data)




