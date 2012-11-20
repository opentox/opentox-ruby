# pending: package dir hack ---------
# CONFIG[:base_dir] = "/home/<user>/opentox-ruby/www"
# PACKAGE_DIR = "/home/<user>/opentox-ruby/r-packages"
package_dir = CONFIG[:base_dir].split("/")
package_dir[-1] = "r-packages"
package_dir = package_dir.join("/")
PACKAGE_DIR = package_dir

require "tempfile"

class Array
  
  def check_uniq
    hash = {}
    self.each do |x|
      raise "duplicate #{x}" if hash[x]
      hash[x] = true
    end
  end
  
end

module OpenTox
  
  class RUtil
    
    @@feats = {}
      
    def initialize
      @r = RinRuby.new(true,false) unless defined?(@r) and @r
      @r.eval ".libPaths('#{PACKAGE_DIR}')"
      @r_packages = @r.pull "installed.packages()[,1]"
      ["sampling","gam","vegan"].each{|l| install_package(l)} #"caret", "smacof", "TunePareto"
      @r.eval "source('#{File.join(Gem.loaded_specs['opentox-ruby'].full_gem_path,'lib/stratification.R')}')"
    end
    
    def quit_r
      begin
        @r.quit
        @r = nil
      rescue
      end
    end
    
    def r
      @r
    end
    
    def package_installed?( package )
      @r_packages.include?(package) 
    end
    
    def install_package( package )
      unless package_installed?(package)
        LOGGER.debug "r-util> installing r-package #{package} to #{PACKAGE_DIR}"
        @r.eval "install.packages('#{package}', repos='http://cran.r-project.org', lib='#{PACKAGE_DIR}')"
      end
    end
    
    # <0 -> array1 << array2
    # 0  -> no significant difference
    # >0 -> array2 >> array1
    def paired_ttest(array1, array2, significance_level=0.95)
      @r.assign "v1",array1
      @r.assign "v2",array2
      @r.eval "ttest = t.test(as.numeric(v1),as.numeric(v2),paired=T)"
      t = @r.pull "ttest$statistic"
      p = @r.pull "ttest$p.value"
      if (1-significance_level > p)
        t
      else
        0
      end
    end
    
    # example: 
    # files = ["/tmp/box.svg","/tmp/box.png"]
    # data = [ [ :method, [4,4,5,5,4,3,2] ], [ :method2, [1,2,3,4,5,4,6] ], [ :asdf, [9,1,8,0,7,1,6] ] ]
    # boxplot(files, data, "comparison1" )
    #
    def boxplot(files, data, title="")
      LOGGER.debug("r-util> create boxplot")
      assign_dataframe("boxdata",data.collect{|e| e[1]}.transpose,nil,data.collect{|e| e[0].to_s})
      plot_to_files(files) do |file|
        @r.eval "boxplot(boxdata,main='#{title}',col=rep(2:#{data.size+1}))"
      end
    end

    # embedds feature values of two datasets into 2D and plots it
    #        
    def feature_value_plot(files, dataset_uri1, dataset_uri2, dataset_name1, dataset_name2,
        features=nil, subjectid=nil, waiting_task=nil)
        
      LOGGER.debug("r-util> create feature value plot")
      d1 = OpenTox::Dataset.find(dataset_uri1,subjectid)
      d2 = OpenTox::Dataset.find(dataset_uri2,subjectid)
      if features
        [d1, d2].each{|d| features.each{|f| raise "feature not included" unless d.features.keys.include?(f)}} 
      else
        raise "different\n#{d1.features.keys.sort.to_yaml}\n#{d2.features.keys.sort.to_yaml}" if 
          (d1.features.keys.sort != d2.features.keys.sort)
        features = d1.features.keys
      end
      raise "at least two features needed" if d1.features.keys.size<2
      waiting_task.progress(25) if waiting_task
      
      df1 = dataset_to_dataframe(d1,0,subjectid,features)
      df2 = dataset_to_dataframe(d2,0,subjectid,features)
      waiting_task.progress(50) if waiting_task
      
      @r.eval "df <- rbind(#{df1},#{df2})"
      @r.eval "split <- c(rep(0,nrow(#{df1})),rep(1,nrow(#{df2})))"
      @r.names = [dataset_name1, dataset_name2]
      LOGGER.debug("r-util> - convert data to 2d")
      #@r.eval "save.image(\"/tmp/image.R\")"
      @r.eval "df.2d <- plot_pre_process(df, method='sammon')"
      waiting_task.progress(75) if waiting_task
      
      LOGGER.debug("r-util> - plot data")
      plot_to_files(files) do |file|
        @r.eval "plot_split( df.2d, split, names, main='Sammon embedding of #{features.size} features',xlab='x',ylab='y')"
      end
    end
    
    # plots a double histogram
    # data1 and data2 are arrays with values, either numerical or categorial (string values)
    # is_numerical, boolean flag indicating value types
    # log (only for numerical), plot logarithm of values
    def double_hist_plot(files, data1, data2, is_numerical, log=false, name1="first", name2="second", title="title", xaxis="x-values")
      LOGGER.debug("r-util> create double hist plot")
      all = data1 + data2
      if (is_numerical)
        @r.eval "double_plot <- function(data1, data2, log=FALSE, names=c('data1','data2'), title='title', xlab='x-values')
        {
          if (log)
          {
            data1 <- log(data1)
            data2 <- log(data2)
            xlab = paste('logarithm of',xlab,sep=' ')
          }
          xlims <- round(c(min(c(min(data1),min(data2))),max(c(max(data1),max(data2)))))
          h <- hist(rbind(data1,data2),plot=F)
          h1 <- hist(data1,plot=F,breaks=h$breaks)
          h2 <- hist(data2,plot=F,breaks=h$breaks)
          xlims = c(min(h$breaks),max(h$breaks))
          ylims = c(0,max(h1$counts,h2$counts))
          xaxps = c(min(h$breaks),max(h$breaks),(length(h$breaks)-1))
          plot(h1, col=rgb(1,0,0,2/4), xlim=xlims, xaxp=xaxps, ylim=ylims,
            main=title, xlab=xlab, ylab='counts' )
          plot(h2, col=rgb(0,1,0,2/4), add=T )
          legend('topleft',names,lty=c(1,1),col=c('red','green'))
        }" 
        @r.assign("data1",data1)
        @r.assign("data2",data2)
        @r.legend = [name1, name2]
      else
        raise "log not valid for categorial" if log
        vals = all.uniq.sort!
        counts1 = vals.collect{|e| data1.count(e)}
        counts2 = vals.collect{|e| data2.count(e)}
        @r.data1 = counts1
        @r.data2 = counts2
        @r.value_names = [name1, name2]
        @r.legend = vals
        @r.eval("data <- cbind(data1,data2)")
      end
      
      plot_to_files(files) do |file|
        if (is_numerical)
          @r.eval "double_plot(data1,data2,log=#{log ? "T":"F"},names=legend,title='#{title}',xlab='#{xaxis}')"
        else
          @r.eval("bp <- barplot(data, beside=T, names.arg=value_names, 
            main='#{title}', col=sort(rep(2:3,length(legend))))") #legend.text=c(legend),
          @r.eval "text(bp, 0, round(data, 1),cex=1,pos=3)"
        end
      end
    end
    
    # stratified splits a dataset into two dataset according to the feature values
    # all features are taken into account unless <split_features> is given
    # returns two datases
    def stratified_split( dataset, metadata={}, missing_values="NA", pct=0.3, subjectid=nil, seed=42, split_features=nil )
      stratified_split_internal( dataset, metadata, missing_values, nil, pct, subjectid, seed, split_features )
    end
    
    # stratified splits a dataset into k datasets according the feature values
    # all features are taken into account unless <split_features> is given
    # returns two arrays of datasets
    def stratified_k_fold_split( dataset, metadata={}, missing_values="NA", num_folds=10, subjectid=nil, seed=42, split_features=nil )
      stratified_split_internal( dataset, metadata, missing_values, num_folds, nil, subjectid, seed, split_features )
    end    
    
    private
    def stratified_split_internal( dataset, metadata={}, missing_values="NA", num_folds=nil, pct=nil, subjectid=nil, seed=42, split_features=nil )
      raise "internal error" if num_folds!=nil and pct!=nil
      k_fold_split = num_folds!=nil
      if k_fold_split
        raise "num_folds not a fixnum: #{num_folds}" unless num_folds.is_a?(Fixnum)
      else
        raise "pct is not a numeric: #{pct}" unless pct.is_a?(Numeric)
      end
      raise "not a loaded ot-dataset (#{dataset.class})" unless dataset.is_a?(OpenTox::Dataset) and dataset.compounds.size>0 and dataset.features.size>0
      raise "missing_values=#{missing_values}" unless missing_values.is_a?(String) or missing_values==0
      raise "subjectid=#{subjectid}" unless subjectid==nil or subjectid.is_a?(String)          
      LOGGER.debug("r-util> apply stratified split to #{dataset.uri}")
      
      df = dataset_to_dataframe( dataset, missing_values, subjectid)
      @r.eval "set.seed(#{seed})"
      str_split_features = ""
      if split_features
        @r.split_features = split_features if split_features
        str_split_features = "colnames=split_features"
      end
      #@r.eval "save.image(\"/tmp/image.R\")"
      
      if k_fold_split
        @r.eval "split <- stratified_k_fold_split(#{df}, num_folds=#{num_folds}, #{str_split_features})"
        split = @r.pull 'split'
        train = []
        test = []
        num_folds.times do |f|
          datasetname = 'dataset fold '+(f+1).to_s+' of '+num_folds.to_s           
          metadata[DC.title] = "training "+datasetname 
          train << split_to_dataset( df, split, metadata, subjectid ){ |i| i!=(f+1) }
          metadata[DC.title] = "test "+datasetname
          test << split_to_dataset( df, split, metadata, subjectid ){ |i| i==(f+1) }
        end
        return train, test
      else
        puts "split <- stratified_split(#{df}, ratio=#{pct}, #{str_split_features})"
        @r.eval "split <- stratified_split(#{df}, ratio=#{pct}, #{str_split_features})"
        split = @r.pull 'split'
        metadata[DC.title] = "Training dataset split of "+dataset.uri
        train = split_to_dataset( df, split, metadata, subjectid ){ |i| i==1 }
        metadata[DC.title] = "Test dataset split of "+dataset.uri
        test = split_to_dataset( df, split, metadata, subjectid ){ |i| i==0 }
        return train, test
      end
    end
    public
    
    # dataset should be loaded completely (use Dataset.find)
    # takes duplicates into account
    # replaces missing values with param <missing_value>
    # returns dataframe-variable-name in R
    def dataset_to_dataframe( dataset, missing_value="NA", subjectid=nil, features=nil )
      LOGGER.debug "r-util> convert dataset to dataframe #{dataset.uri}"
      
      # use either all, or the provided features, sorting is important as col-index := features
      if features
        features.sort!
      else
        features = dataset.features.keys.sort
      end

      # values into 2D array, then to dataframe
      d_values = []
      dataset.compounds.size.times do |c_idx|
        c_values = []
        features.each do |f|
          v = dataset.data_entry_value(c_idx,f)
          v = missing_value if v==nil
          c_values << v
        end
        d_values << c_values
      end  
      df_name = "df_#{dataset.uri.split("/")[-1].split("?")[0]}"
      
      compound_names = dataset.compounds.size.times.collect{|idx| dataset.compounds[idx]+"$"+idx.to_s}
      assign_dataframe(df_name,d_values,compound_names,features)
      
      # set dataframe column types accordingly
      f_count = 1 #R starts at 1
      features.each do |f|
        feat = OpenTox::Feature.find(f,subjectid)
        nominal = feat.metadata[RDF.type].to_a.flatten.include?(OT.NominalFeature)
        if nominal
          @r.eval "#{df_name}[,#{f_count}] <- as.character(#{df_name}[,#{f_count}])"
        else
          @r.eval "#{df_name}[,#{f_count}] <- as.numeric(#{df_name}[,#{f_count}])"
        end
        f_count += 1
      end
      #@r.eval "head(#{df_name})"
      
      # store compounds, and features (including metainformation)
      @@feats[df_name] = {}
      features.each do |f|
        @@feats[df_name][f] = dataset.features[f]
      end
      df_name
    end
    
    # converts a dataframe into a dataset (a new dataset is created at the dataset webservice)
    # this is only possible if a superset of the dataframe was created by dataset_to_dataframe (metadata and URIs!)
    def dataframe_to_dataset( df, metadata={}, subjectid=nil )
      dataframe_to_dataset_indices( df, metadata, subjectid, nil)
    end
    
    private
    def dataframe_to_dataset_indices( df, metadata={}, subjectid=nil, compound_indices=nil )
      raise unless @@feats[df].size>0
      values, compound_names, features = pull_dataframe(df)
      compounds = compound_names.collect{|c| c.split("$")[0]}
      features.each{|f| raise unless @@feats[df][f]}
      dataset = OpenTox::Dataset.create(CONFIG[:services]["opentox-dataset"],subjectid)
      dataset.add_metadata(metadata)
      LOGGER.debug "r-util> convert dataframe to dataset #{dataset.uri}"
      features.each{|f| dataset.add_feature(f,@@feats[df][f])}
      feature_numeric = features.size.times.collect do |c|
        feat = OpenTox::Feature.find(features[c],subjectid)
        feat.metadata[RDF.type].to_a.flatten.include?(OT.NumericFeature) 
      end
      compounds.size.times do |r|
        if compound_indices==nil or compound_indices.include?(r)
          dataset.add_compound(compounds[r])          
          features.size.times do |c|
            dataset.add_data_entry(compounds[r],features[c],feature_numeric[c] ? values[r][c].to_f : values[r][c]) if values[r][c]!="NA"
          end
        end
      end
      dataset.save(subjectid)
      dataset
    end    
    
    def split_to_dataset( df, split, metadata={}, subjectid=nil )
      indices = []
      split.size.times{|i| indices<<i if yield(split[i]) }
      dataset = dataframe_to_dataset_indices( df, metadata, subjectid, indices )
      LOGGER.debug("r-util> split into #{dataset.uri}, c:#{dataset.compounds.size}, f:#{dataset.features.size}")
      dataset
    end
    
    def pull_dataframe(df)
      tmp = File.join(Dir.tmpdir,Time.new.to_f.to_s+"_"+rand(10000).to_s+".csv")
      @r.eval "write.table(#{df},file='#{tmp}',sep='#')"
      res = []; compounds = []; features = []
      first = true
      file = File.new(tmp, 'r')
      file.each_line("\n") do |row|
        if first
           features = row.chomp.split("#").collect{|e| e.gsub("\"","")}
           first = false
        else
           vals = row.chomp.split("#").collect{|e| e.gsub("\"","")}
           compounds << vals[0]
           res << vals[1..-1]
        end
      end
      begin File.delete(tmp); rescue; end
      return res, compounds, features
    end
    
    def assign_dataframe(df,input,rownames,colnames)
      rownames.check_uniq if rownames
      colnames.check_uniq if colnames
      tmp = File.join(Dir.tmpdir,Time.new.to_f.to_s+"_"+rand(10000).to_s+".csv")
      file = File.new(tmp, 'w')
      input.each{|i| file.puts(i.collect{|e| "\"#{e}\""}.join("#")+"\n")}  
      file.flush
      @r.rownames = rownames if rownames
      @r.colnames = colnames
      @r.eval "#{df} <- read.table(file='#{tmp}',sep='#',"+
        "#{rownames ? "row.names=rownames" : ""},col.names=colnames,check.names=F)"
      begin File.delete(tmp); rescue; end
    end
    
    def plot_to_files(files)
      files.each do |file|
        if file=~/(?i)\.svg/
          @r.eval("svg('#{file}',10,8)")
        elsif file=~/(?i)\.png/
          @r.eval("png('#{file}')")
        else
          raise "invalid format: "+file.to_s
        end
        yield file
        LOGGER.debug "r-util> plotted to #{file}"
        @r.eval("dev.off()")
      end
    end
  end
end


