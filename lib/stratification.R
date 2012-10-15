
# 	bbrc dataset 249 x 530
#load("/home/martin/workspace/ValidationExperiments/bbrc_image.R")
#data=df_12117 

# 250 compounds, 2 features
#load("/home/martin/workspace/ValidationExperiments/strat_pics/image.R")
#data=df_11306
#data=cbind(data[1],data[3])

# 1000 compounds, 226 features
#load("/home/martin/tmp/image_12171.R")
#data=df_12171

is_really_numeric <- function(data)
{
	for (i in 1:ncol(data))
	{
		if (!is.numeric(data[,i]))
			return(FALSE)
	}
	TRUE
}

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
  #save.image("/tmp/to_process.R")
	
  rownames(data) <- NULL
  data.num <- as.data.frame(data)
  if (!is.null(colnames))
  {
    data.num = subset(data.num, select = colnames)
	## if (ncol(data.num)==0)
	## {
	##     print("ERROR no columns left after selecting colnames")
	##     print("colnames before selection:")
	##     print(colnames(data))
	##     print("colnames to be selected:")
	##     print(colnames)
	##     stop("no columns left after selecting colnames")
	## }
  }
  if (!is_really_numeric(data.num))
    data.num = nominal_to_binary(data.num)
  if(any(is.na(data.num)))
  {
	print("WARNING data contains NAs, if you convert to dataframe you can specify the missing values value.")
	suppressPackageStartupMessages(require("gam"))
   	data.repl = na.gam.replace(data.num)
  }
  else
  	data.repl = data.num
  data.without_constant_vals = remove_var_zero(data.repl)
  if(ncol(data.repl)>ncol(data.without_constant_vals))
	  print(paste("removed ",(ncol(data.repl)-ncol(data.without_constant_vals)),"/",ncol(data.repl)," features because of no variance"))
  if (ncol(data.without_constant_vals)==0)
	  stop("ERROR no columns left after processing")
  data.without_constant_vals
}

# depends on random seed
cluster_knn <- function( data, min=10, max=15) #, method="kmeans" )
{
  suppressPackageStartupMessages(require("vegan"))
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

dynamic_dist <- function(data, feature_type)
{
	binary = check_binary(data, feature_type)
	if (!binary)
	{
		print(paste("distance used: Euclidean"))
		d <- dist(data, method="euclidean")
	}
	else
	{
		print(paste("distance used: Jaccard (Tanimoto)"))
		suppressPackageStartupMessages(require("proxy"))
		d <- dist(data, method="Jaccard")
	}
	m <- as.matrix(d) 
	diag(m)=0 # do this manually, as it is set to NA by the Tanimoto method
	list(dist=d,matrix=m)
}

#deterministic!
cluster_hierarchical <- function(data, feature_type, hclust_method="ward", max_num_clusters=15, deep_split=4, ...)
{
	max_num_clusters <- min(max_num_clusters,nrow(unique(data)))
	max_num_clusters <- min(max_num_clusters,nrow(data)-1)
	min_size <- round_it(nrow(data)/max_num_clusters)
	d <- dynamic_dist(data, feature_type)	
	suppressPackageStartupMessages(require("dynamicTreeCut"))
	fit <- hclust(d$dist, method=hclust_method)
	cut = cutreeDynamic(fit, method = "hybrid", distM = d$matrix, minClusterSize = min_size, deepSplit=deep_split, ...)
	print(paste("dynamicTreeCut clustering, minClusterSize:",min_size," result:"))
	print(table(cut))
	#print(cut)
	#as.data.frame(cut)
	cut
}

# input: classes [ 1 2 4 3 1 2 4 3 2 1 2 3 3 2 ]
# input : ratio 0.1
# output: split  [ 1 0 0 1 0 0 0 0 0 1 0 0 0 0 ]
explicit_stratified_split <- function( class, ratio=0.3 )
{
#	print(class)
	num_instances = length(class)
#	print(num_instances)

	min_class=min(class)
	max_class=max(class)
	
	# create a list with indexes of each class
	class_idx <- list()
	for (i in min_class:max_class)
		class_idx[[i+1]] <- vector() 
	for (i in 1:num_instances)
		class_idx[[class[i]+1]] <- c(class_idx[[class[i]+1]],i)
	#print(class_idx)
	# shuffle each index list
	for (i in min_class:max_class)
		class_idx[[i+1]] <- sample(class_idx[[i+1]]) 
 #   print(class_idx)
	
	selected <- vector()
	sum_total <- 0
	
	classes = sample(min_class:max_class)
	# iterate of classes in random order
	for (i in classes)
	{
		# count the numbers of selections that should be made
		sum_total <- sum_total + length(class_idx[[i+1]])
		selected_total <- round_it(sum_total*ratio)
	#	print(selected_total)
		# count the number of new selections
		selected_new <- selected_total - length(selected)
	#	print(selected_new)
		# add selections
		selected <- c(selected, class_idx[[i+1]][1:selected_new])
	#	print(selected)
	#	print("")
	}
	
	# convert selected indexs to split array
	split <- array(0,num_instances)
	for(i in selected)
		split[i] <- 1
	as.vector(split)
}

random_split <- function( data, ratio=0.3 )
{
	num <- round_it(nrow(data)*ratio)
	selected <- sample(rep(1:nrow(data)))[1:num]
	# convert selected indexs to split array
	split <- array(0,nrow(data))
	for(i in selected)
		split[i] <- 1
	as.vector(split)
}

## is_binary <- function( data )
## {
##     #print(length(data[data==0]))
##     #print(length(data[data==1]))
##     #print(length(data[data==0])+length(data[data==1]))
##     #print(nrow(data))
##     #print(ncol(data))
##     #print(nrow(data)*ncol(data))
##     #print(head(data))
##     length(data[data==0])+length(data[data==1]) == nrow(data)*ncol(data)
## }

check_binary <- function( data, feature_type)
{
	binary <- length(data[data==0])+length(data[data==1]) == nrow(data)*ncol(data) #=is_binary(data)
	if (feature_type=="numerical")
	{
		if (binary)
			stop("data is binary, but specified feature_type is numerical")
	}
	else if (feature_type=="binary")
	{
		if (!binary)
			stop("data is not binary, but specified feature_type is binary")
	}
	else
		stop(paste("unknown feature type",feature_type))
	binary
}


contra_stratified_split <- function( data, feature_type, ratio=0.3, colnames=NULL ) #, samplesize=10 )
{
	data.processed = as.matrix(process_data( data, colnames ))
	print(paste("contra strat split, #features:",ncol(data.processed),"/",ncol(data),", ratio:",ratio))

	binary=check_binary(data.processed, feature_type)
	
	##save.image("/tmp/contra_splitting.R")
	
	
	print(paste("data is binary: ",binary))
	
	if (!binary)
	{
		data.processed = as.matrix(pca_reduce_features(data.processed))
		print(paste("#features reduced with pca: ",ncol(data.processed)))
	}
	else
	{
		#do nothing
	}
	
	samplesize = nrow(data.processed)/20
	print(paste("samplesize",samplesize))
	if (nrow(data.processed)<=samplesize)
	   sample = as.vector(1:nrow(data.processed))
	else
	   sample = random_split(data,ratio=samplesize/nrow(data))
	sample_data = data.frame()
	orig_idx = c()
	for(i in 1:nrow(data.processed))
	    if(sample[i]==1)
	    {
	        sample_data = rbind(sample_data, data.processed[i,])
	        orig_idx = cbind(orig_idx,i)
	    }
	#print("head(sample data) (random sampled with fixed samplesize)")
	#print(head(sample_data))
	print("orig_idx (orinal indices of sample data in orig data)")
	print(as.vector(orig_idx))
	
	m <- dynamic_dist(sample_data, feature_type)$matrix
	
	max_dist = 0
	max_dist_idx = -1
	for(i in 1:nrow(sample_data))
	{
	    m_dist <- sum(m[,i])
	    if (max_dist_idx==-1 || m_dist>max_dist)
	    {
	        max_dist = m_dist
	        max_dist_idx = i
	    }
	}
	print("max-dist (in sample distance matrix)")
	print(max_dist)
	print("max-dist-id")
	print(max_dist_idx)
	
	anti_strat_center = orig_idx[max_dist_idx]
	print("anti strat center (orignal index of max-distindex)")
	print(anti_strat_center)
	
	#anti_strat_center = sample(1:nrow(data.processed),1)
	
	dist = array(0,nrow(data.processed))
	if (!binary)
	{
		for(i in 1:nrow(data.processed))
			dist[i] = sqrt(sum((data.processed[i,] - data.processed[anti_strat_center,]) ^ 2))
	}
	else
	{
		for(i in 1:nrow(data.processed))
			dist[i] = as.vector(dist(rbind(data.processed[i,],data.processed[anti_strat_center,]),method="Jaccard"))
	}
	
	print("head(dist) (distance from all data to anti strat center)")
	print(head(dist))
	max = max(dist)
	print("max (maxium distance)")
	print(max)
	
	not_nil = 0
	prob = array(0,nrow(data.processed))
	#raw = array(0,nrow(data.processed))
	for(i in 1:nrow(data.processed))
	{
		prob[i] = (1-dist[i]/max) ^ 100		
		#raw[i] = 1-dist[i]/max
	    if(prob[i]!=0)
			not_nil = not_nil+1
	}
	
	
	print("head(prop) (convert distance into propability)")
	print(head(prob))
	sum_prop = sum(prob)
	print("sum prop (sum of propabilities)")
	print(sum_prop)
	
	num_sel = round_it(nrow(data.processed)*ratio)
	print(paste("selecting",num_sel,"from",nrow(data.processed),"according to prob ( prob!=0 for",not_nil,")"))
	if (num_sel>not_nil)
	{
	    for(i in 1:nrow(data.processed))
			if(prob[i]==0)
				prob[i]=.Machine$double.xmin
	}	
	
	selected = 	sample(1:nrow(data.processed),num_sel,replace=FALSE,prob=prob)
	split <- array(0,nrow(data.processed))
	for(i in selected)
	{
		#print(paste("sel ",i,"  raw ",raw[i]," dist ",prop[i]))
		split[i] <- 1
	}
	split = as.vector(split)
	
	cl <- array(0,nrow(data.processed))
	cl[anti_strat_center] = 1
	
	list(split=split,cluster=cl)
}

stratified_split <- function( data, feature_type, ratio=0.3, method="cluster_knn", method_2="samplecube", colnames=NULL ) #, preprocess="none"
{
	data.processed = as.matrix(process_data( data, colnames ))
	print(paste("strat split, method: ",method," #features: ",ncol(data.processed)," ratio: ",ratio)) #," preprocess: ",preprocess))
	cl = NULL
	
	binary = check_binary(data.processed, feature_type)
	
	if (!binary)
	{
		data.processed = as.matrix(pca_reduce_features(data.processed))
		print(paste("#features reduced with pca: ",ncol(data.processed)))
	}
	
    if (method == "samplecube")
    {
	  suppressPackageStartupMessages(require("sampling"))
      # adjust ratio to make samplecube return exact number of samples
      ratio = round_it(nrow(data.processed)*ratio)/nrow(data.processed)
      pik = rep(ratio,times=nrow(data.processed))
      data.strat = cbind(pik,data.processed)
      split = samplecube(data.strat,pik,order=2,comment=F)
    }
    else if (method == "cluster_knn")
    {
      cl = as.vector(cluster_knn(data.processed))
#      suppressPackageStartupMessages(require("caret"))
#      res = createDataPartition(cl,p=ratio)
#      split = rep(1, times=nrow(data))
#      for (j in 1:nrow(data))
#        if ( is.na(match(j,res$Resample1)) )
#          split[j]=0
#      split
	   
       split = stratified_split(cl,ratio,"samplecube")$split
    }
	else if (method == "cluster_hierarchical") 
	{
		cl = as.vector(cluster_hierarchical(data.processed, feature_type))
		
	    #suppressPackageStartupMessages(require("caret"))
        #res = createDataPartition(cl,p=ratio)
        #split = rep(1, times=nrow(data))
        #for (j in 1:nrow(data))
        #if ( is.na(match(j,res$Resample1)) )
        #  split[j]=0
		
	   if (method_2 == "samplecube")
	     split = stratified_split(cl,ratio,"samplecube")$split
	   else if (method_2 == "explicit")
	     split = explicit_stratified_split(cl,ratio)
	   else
		   stop("unknown method")
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
	##     suppressPackageStartupMessages(require("sampling"))
	##     stratified_split(cl,ratio,"samplecube")
	##   }
	## }
    else
      stop("unknown method")
  
    list(split=split, cluster=cl)
}

anti_stratified_split_OLD <- function( data, ratio=0.3, colnames=NULL)
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
  cl = cluster_knn(data.processed, num_c, num_c)
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
  list(split=as.vector(split),cluster=cl)
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
	  suppressPackageStartupMessages(require("sampling"))
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
	suppressPackageStartupMessages(require("TunePareto"))
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

sammon_duplicates <- function( data, feature_type, ... ) { 
  di <- duplicate_indices(data)
  #print("duplicate indices")
  #print(di)
  u <- unique(data) 
  print(paste('unique data points',nrow(u),'of',nrow(data)))
  if(nrow(u) <= 4) stop("number of unqiue datapoints <= 4")
  distance <- dynamic_dist(u, feature_type)$dist
  points_unique <- sammon(distance, ...)$points
  #print("points unique")
  #print(points_unique)
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
  
  #print("pca reduce features2")
  
  data.var <- data.pca$sdev^2
  
  #print("pca reduce features3")
  
  explained_variance <- cumsum(data.var)/sum(data.var)
  
  print('explained variance:')
  print(explained_variance)
  for (i in 1:length(explained_variance))
	  if (i>=max_num_features || explained_variance[i]>=min_explained_variance)
		  break
  as.data.frame(data.pca$x)[1:i]
}

plot_pre_process <- function( data, feature_type, method="pca" )
{
  data.processed = process_data( data )
  if (method == "pca")
  {
    data.pca <- prcomp(data.processed, scale=TRUE)
    as.data.frame(data.pca$x)[1:2]
  }
  else if (method == "smacof")
  {
	suppressPackageStartupMessages(require("smacof"))
    data.emb <- smacofSym(dist(data.processed, method = "euclidean"), ndim=2, verbose=T)
    data.emb$conf
  }
  else if (method == "sammon")
  {
	suppressPackageStartupMessages(require("MASS"))
    sammon_duplicates(data.processed, feature_type, k=2)
  }
  else
    stop("unknown method")
}


plot_split <- function( data, color_idx=NULL, circle_idx=NULL, ... )
{ 
	if (ncol(data)!=2 || !is.numeric(data[,1]) || !is.numeric(data[,2]))
		stop("data not suitable for plotting, plot_pre_process() first")
	
	plot( NULL, xlim = extendrange(data[,1]), ylim = extendrange(data[,2]), ... )
	if (is.null(names))
		names <- c("split 1","split 2")
	#colos = as.double(rep(2:(max(color_idx)+2)))
	#legend("topleft",names,pch=2,col=colos)
	
	#legend("topleft",names[1],pch=as.double(c(1)),col="red")
	
	## for (j in max(color_idx):0)
	## {
	##   for (k in max(shape_idx):0)
	##   {
	##     set = c()
	##     for (i in 1:nrow(data))
	##       if (color_idx[i]==j && shape_idx[i]==k)
	##         set = c(set,i)
	##     points(data[set,], pch = 15+k, col=(j+2))
	##   }
	## }
	
	
	col_offset = 2
	if(is.null(color_idx))
	{
		color_idx <- array(0,nrow(data))
		col_offset=3
	}
	for (j in 0:max(color_idx))
	{
		set = c()
		for (i in 1:nrow(data))
			if (color_idx[i]==j)
				set = c(set,i)
		points(data[set,], pch = 19, cex=1, col=(max(color_idx)-j)+col_offset)
	}
	if (!is.null(circle_idx))
	{
		set = c()
		for (i in 1:nrow(data))
			if (circle_idx[i]==1)
				set = c(set,i)#
		
		#print(set)
		#print(data[set,])
		
		points(matrix(data[set,],ncol=2), pch = 1, cex=1, col=1)
		points(matrix(data[set,],ncol=2), pch = 1, cex=1.8, col=1)
		#points(data[set,], pch = 1, cex=1, col=1)
		#points(data[set,], pch = 1, cex=1.8, col=1)
		
		## for (j in max(color_idx):0)
		## {
		##     set = c()
		##     for (i in 1:nrow(data))
		##         if (color_idx[i]==j && circle_idx[i]==1)
		##             set = c(set,i)
		##     points(data[set,], pch = 19, cex=1, col=1)
		##     points(data[set,], pch = 20, cex=1, col=(max(color_idx)-j)+col_offset)
		##     points(data[set,], pch = 1, cex=2, col=1)
		## }
		
	}
}

superboxplot  <- function( data, ..., significance_level=0.95, paired=T, test_closer_to_zero=F, alg_info=T )
{
	b <- boxplot(data,...) #,col=rep(2:(ncol(data)+1)))
	
	if (alg_info)
	{
		#print mean and stdev
		for (i in 1:ncol(data))
		{
			med <- sprintf("%.3f",b$stats[,i][3])
			stdev <- sprintf("%.3f",sqrt(var(data[,i])))
			mtext(paste(med,"+-",stdev),side=1,at=i,line=2)
		}
		
		#print significance tests
		if (nrow(data)>10 && significance_level>0 && significance_level<=1)
		{
			sig <- array(list(),c(ncol(data),ncol(data)))
			sig_var <- array(list(),c(ncol(data),ncol(data)))
			for (i in 1:(ncol(data)-1))
			{
				for (j in (i+1):ncol(data))
				{
					if (test_closer_to_zero==F)
						ttest = t.test(data[,i],data[,j],paired=paired)
					else
						ttest = ttest_closer_to_zero(data[,i],data[,j],paired=paired)
					if ( !is.na(ttest$p.value) && 1-significance_level > ttest$p.value)
					{
						sig[i,j] = T
						sig[j,i] = T
					}
					else
					{
						sig[i,j] = F
						sig[j,i] = F
					}
					
					ftest = var.test(data[,i],data[,j])
					if ( !is.na(ftest$p.value) && 1-significance_level > ftest$p.value)
					{
						sig_var[i,j] = T
						sig_var[j,i] = T
					}
					else
					{
						sig_var[i,j] = F
						sig_var[j,i] = F
					}
				}
			}
			for (i in 1:ncol(data))
			{
				## s <- ""
				## for (j in 1:ncol(data))
				## {
				##     if (i == j)
				##         s <- paste( s, "-" ,sep="")
				##     else if (sig[i,j]==T)
				##         s <- paste( s, "X" ,sep="")
				##     else
				##         s <- paste( s, "0" ,sep="")
				## }
				
				s<-""
				bigger <- ""
				for (j in 1:ncol(data))
				{
					#print(paste(i,j))
					if (i!=j && sig[i,j]==T)
					{
						if (test_closer_to_zero==F && b$stats[,i][3] > b$stats[,j][3]) 
							bigger <- paste(bigger,j,sep=",")
						else if (test_closer_to_zero==T && abs(b$stats[,i][3]) > abs(b$stats[,j][3])) 
							bigger <- paste(bigger,j,sep=",")
					}
				}
				if (nchar(bigger)>0)
				{
					s<-"med"
					bigger <- substring(bigger, 2)
					bigger <- paste(">(",bigger,")",sep="")
					s <- paste(s,bigger,sep=" ")
				}
				smaller <- ""
				for (j in 1:ncol(data))
				{
					#print(paste(i,j))
					if (i!=j && sig[i,j]==T)
					{
						if (test_closer_to_zero==F && b$stats[,i][3] < b$stats[,j][3]) 
							smaller <- paste(smaller,j,sep=",")
						else if (test_closer_to_zero==T && abs(b$stats[,i][3]) < abs(b$stats[,j][3])) 
							smaller <- paste(smaller,j,sep=",")
					}
				}
				if (nchar(smaller)>0)
				{
					if(nchar(bigger)==0)
						s<-"med"
					smaller <- substring(smaller, 2)
					smaller <- paste("<(",smaller,")",sep="")
					s <- paste(s,smaller,sep=" ")
				}
				
				bigger <- ""
				for (j in 1:ncol(data))
				{
					if (i!=j && sig_var[i,j]==T && var(data[,i]) > var(data[,j])) 
						bigger <- paste(bigger,j,sep=",")
				}
				if (nchar(bigger)>0)
				{
					s <- paste(s,"var",sep=" ")
					bigger <- substring(bigger, 2)
					bigger <- paste(">(",bigger,")",sep="")
					s <- paste(s,bigger,sep=" ")
				}
				smaller <- ""
				for (j in 1:ncol(data))
				{
					#print(paste(i,j))
					if (i!=j && sig_var[i,j]==T && var(data[,i]) < var(data[,j])) 
						smaller <- paste(smaller,j,sep=",")
				}
				if (nchar(smaller)>0)
				{
					if(nchar(bigger)==0)
						s <- paste(s,"var",sep=" ")
					smaller <- substring(smaller, 2)
					smaller <- paste("<(",smaller,")",sep="")
					s <- paste(s,smaller,sep=" ")
				}			
				
				
				mtext(s,side=1,at=i,line=3)
			}
		}
	}
	
	#print(sig)
}

ttest_closer_to_zero <- function( x, y, ... )
{
	prep = pre_process_ttest_closer_to_zero(x,y)
	t.test(prep$x, prep$y, ...)	
}

pre_process_ttest_closer_to_zero  <- function( x, y )
{
	medX = median(x)
	medY = median(y)
	if((medX > 0 && medY > 0)|| (medX < 0 && medY < 0))
	{
		list(x=x,y=y)
	}
	else
	{
		pos_data = x
		neg_data = y
		if (medX<medY)
		{
			pos_data = y
			neg_data = x
		}
		shift = min(abs(medX),abs(medY))
		pos_data = pos_data - shift
		neg_data = neg_data + shift
		if (medX>medY)
			list(x=pos_data,y=neg_data)
		else
			list(x=neg_data,y=pos_data)
	}
}

split_plot <- function(ratio=0.33)
{
	print("splitting")
	split = contra_stratified_split(data, ratio=ratio)
	#print(split$cluster)
	print("plotting")
	plot_split(plot_data,circle_idx=split$cluster,color_idx=split$split)
}

norm_test <- function(data.orig)
{
	print(paste("cols: ",ncol(data.orig)))
	data <- remove_var_zero(data.orig)
	print(paste("cols var != zero: ",ncol(data)))	
	for (j in 1:ncol(data))
	{
		col = data[,j]
		if(length(col)>5000)
			col = sample(col,5000)
		result = tryCatch({
			w = shapiro.test(col)$statistic
			w_log = shapiro.test(log(col))$statistic
			if (w<0.60 && w_log>0.9)
			{
				print(paste("replace with log","new-index:",j,"old-index:",which(colnames(data.orig)==colnames(data)[j]),colnames(data)[j],w,w_log))
				colnames(data)[j] <- paste("log",colnames(data)[j],sep="_")
				data[,j] <- log(data[,j])
			}
		}, error = function(err){
			#print(paste(j,colnames(d)[j],err))
		}, warning = function(warn){
			#print(paste("warning",warn))
		})
	}
	data
}


#a<-matrix(rnorm(100, mean=50,  sd=4), ncol=5)
#b<-matrix(rnorm(5000, mean=0, sd=30), ncol=5)
#data<-b
#data<-rbind(a,b)
#c<-matrix(rnorm(50, mean=-50, sd=2), ncol=5)
#data<-rbind(data,c)
#data=iris[1:4]


#data = data[1:10,]

#require("datasets")
#data=ChickWeight
#data=freeny
#data=sunspots
#data=faithful
#pca_reduce_features(data)
#data=ethanol
#split = stratified_k_fold_split(data, num_folds=3)

#set.seed(as.numeric(Sys.time()))
#split = contra_stratified_split(data, ratio=0.1)
#split = stratified_split(data, ratio=0.1, method="cluster_hierarchical")
#split = stratified_split(data, ratio=0.1, method="cluster_hierarchical")

#split = stratified_split(data, ratio=0.1,preprocess = "pca", method="cluster_hierarchical", method_2="explicit")
#split = stratified_split(data, ratio=0.1,preprocess = "pca", method="cluster_hierarchical")

#cluster = cluster(ethanol,3,3)

#print(split)
#print(sum(split)

#print(split$cluster)
#print(which(split$cluster==1))

#plot_data = plot_pre_process(data, method="sammon")

#print(nrow(plot_data))
#print(ncol(plot_data))
#print(head(plot_data))
#print(plot_data[which(split$cluster==1),])


#plot_data = data
#plot_split(plot_data,color_idx=split$split)
#plot_split(plot_data,circle_idx=split$split,color_idx=split$cluster)
#plot_split(plot_data,circle_idx=split$cluster,color_idx=split$split)



#cl = cluster(data)



