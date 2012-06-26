
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

remove_var_zero <- function(data, useNA = 'ifany')
{
	out <- apply(data, 2, function(x) {length(table(x, useNA = useNA))})
	if (any(out==1))
		data[,-which(out==1)]
	else
		data
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
  data.without_constant_vals = remove_var_zero(data.repl)
  data.without_constant_vals
}

# depends on random seed
cluster_knn <- function( data, min=10, max=15) #, method="kmeans" )
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

#deterministic!
cluster_hierarchical <- function(data, hclust_method="ward", max_num_clusters=15, deep_split=4, ...)
{
	max_num_clusters <- min(max_num_clusters,nrow(unique(data)))
	max_num_clusters <- min(max_num_clusters,nrow(data)-1)
	min_size <- round_it(nrow(data)/max_num_clusters)
	
	require("dynamicTreeCut")
	d <- dist(data, method = "euclidean")
	fit <- hclust(d, method=hclust_method)
	cut = cutreeDynamic(fit, method = "hybrid", distM = as.matrix(d), minClusterSize = min_size, deepSplit=deep_split, ...)
	print(paste("dynamicTreeCut clustering, minClusterSize:",min_size," result:"))
	print(table(cut))
	#print(cut)
	#as.data.frame(cut)
	cut
}

stratified_split <- function( data, ratio=0.3, method="cluster_knn", colnames=NULL, preprocess="none" )
{
	data.processed = as.matrix(process_data( data, colnames ))
	print(paste("strat split, method: ",method," #features: ",ncol(data.processed)," ratio: ",ratio," preprocess: ",preprocess))
	
	if (preprocess == "none")
	{
		#do nothing
	}	
	else if (preprocess == "pca")
	{
		data.processed = as.matrix(pca_reduce_features(data.processed))
		print(paste("#features reduced with pca: ",ncol(data.processed)))	
	}
	else
		stop("unknown preprocess")
	
    if (method == "samplecube")
    {
      require("sampling")
      # adjust ratio to make samplecube return exact number of samples
      ratio = round_it(nrow(data.processed)*ratio)/nrow(data.processed)
      pik = rep(ratio,times=nrow(data.processed))
      data.strat = cbind(pik,data.processed)
      samplecube(data.strat,pik,order=2,comment=F)
    }
    else if (method == "cluster_knn")
    {
      cl = cluster_knn(data.processed)
#      require("caret")
#      res = createDataPartition(cl,p=ratio)
#      split = rep(1, times=nrow(data))
#      for (j in 1:nrow(data))
#        if ( is.na(match(j,res$Resample1)) )
#          split[j]=0
#      split
      stratified_split(cl,ratio,"samplecube")
    }
	else if (method == "cluster_hierarchical") 
	{
		cl = cluster_hierarchical(data.processed)
		stratified_split(cl,ratio,"samplecube")
	}
	## else if (method == "cluster2")
	## {
	##   if (ratio > 0.5)
	##   {
	##     num = round_it(nrow(data.processed)*(1-ratio))
	##     swap=TRUE
	##   }
	##   else
	##   {
	##     num = round_it(nrow(data.processed)*ratio)
	##     swap=FALSE
	##   }
	##   cl <- cluster(data.processed, num, num)  
	## 
	##   if(max(cl)==num)
	##   {
	##     split <- array(1:nrow(data))
	##     check <- array(0,num)
	##     count = 0
	##     for(j in sample(array(nrow(data))))
	##     {
	##        if (count<num && check[cl[j]]==0)
	##        {
	##          split[j] = 1
	##          check[cl[j]] = 1
	##          count=count+1
	##        }
	##        else
	##          split[j] = 0
	##     } 
	##     if (swap)
	##       for(j in 1:nrow(data))
	##         split[j] = 1-split[j]
	##     as.vector(split)
	##   }
	##   else
	##   {
	##     require("sampling")
	##     stratified_split(cl,ratio,"samplecube")
	##   }
	## }
    else
      stop("unknown method")
}

anti_stratified_split <- function( data, ratio=0.3, colnames=NULL)
{
  if (ratio > 0.5)
  {
    ratio = 1-ratio
    swap = TRUE
  }
  else
    swap = FALSE
  data.processed = as.matrix(process_data( data, colnames ))
  print(paste("anti-split using #features: ",ncol(data.processed)))
  num_c = floor(1/ratio)
  cl = cluster(data.processed, num_c, num_c)
  #print(cl)
  idx = -1
  min = 1000000
  num = round_it(nrow(data)*ratio)
  for(j in 1:max(cl))
  {
    cl_size = length(subset(cl, cl==j))
    if (cl_size<min && cl_size>=num)
    {
      idx = j
      min = cl_size
    }
  }
  split <- array(1:nrow(data))
  count = 0
  for(j in sample(array(nrow(data))))
  {
     if (count<num && cl[j]==idx)
     {
       split[j] = 1
       count=count+1
     }
     else
       split[j] = 0
     
  } 
  if (swap)
    for(j in 1:nrow(data))
      split[j] = 1-split[j]
  #print(split)
  as.vector(split)
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

pca_reduce_features <- function( data, max_num_features=1000, min_explained_variance=1 )
{
  print("pca reduce features")
  #save.image("/tmp/image.R")

  #data.processed = process_data( data )
  data.pca <- prcomp(data,scale=TRUE)
  
  print("pca reduce features2")
  
  data.var <- data.pca$sdev^2
  
  print("pca reduce features3")
  
  explained_variance <- cumsum(data.var)/sum(data.var)
  
  print('explained variance:')
  print(explained_variance)
  for (i in 1:length(explained_variance))
	  if (i>=max_num_features || explained_variance[i]>=min_explained_variance)
		  break
  as.data.frame(data.pca$x)[1:i]
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

plot_split <- function( data, color_idx, names=NULL, shape_idx=color_idx, ... )
{
  if (ncol(data)!=2 || !is.numeric(data[,1]) || !is.numeric(data[,2]))
    stop("data not suitable for plotting, plot_pre_process() first")

  plot( NULL, xlim = extendrange(data[,1]), ylim = extendrange(data[,2]), ... )
  if (is.null(names))
    names <- c("split 1","split 2")
  colos = as.double(rep(2:(max(split)+2)))
  legend("topleft",names,pch=2,col=colos)

  for (j in max(color_idx):0)
  {
    for (k in max(shape_idx):0)
    {
      set = c()
      for (i in 1:nrow(data))
        if (color_idx[i]==j && shape_idx[i]==k)
          set = c(set,i)
      points(data[set,], pch = 15+k, col=(j+2))
    }
  }
}

#a<-matrix(rnorm(100, mean=50,  sd=4), ncol=5)
#b<-matrix(rnorm(5000, mean=0, sd=10), ncol=5)
#data<-rbind(a,b)
#c<-matrix(rnorm(50, mean=-50, sd=2), ncol=5)
#data<-rbind(data,c)
#data=iris
#pca_reduce_features(data)
#data=ethanol
#split = stratified_k_fold_split(data, num_folds=3)
#split = anti_stratified_split(data, ratio=0.75)

#split = stratified_split(data, ratio=0.9,preprocess = "pca")

#cluster = cluster(ethanol,3,3)

#print(split)
#print(sum(split))
#plot_split(plot_pre_process(data, method="pca"),split,c("training","test"))



#cl = cluster(data)




