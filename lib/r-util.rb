# pending: package dir hack ---------
# CONFIG[:base_dir] = "/home/<user>/opentox-ruby/www"
# PACKAGE_DIR = "/home/<user>/opentox-ruby/r-packages"
package_dir = CONFIG[:base_dir].split("/")
package_dir[-1] = "r-packages"
package_dir = package_dir.join("/")
PACKAGE_DIR = package_dir

require "tempfile"
require "statsample"

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
      ["sampling","gam","vegan","dynamicTreeCut"].each{|l| install_package(l)} #"caret", "smacof", "TunePareto"
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
    def ttest(array1, array2, paired, significance_level=0.95)
      @r.assign "v1",array1
      @r.assign "v2",array2
      raise if paired && array1.size!=array2.size
      @r.eval "ttest = t.test(as.numeric(v1),as.numeric(v2),paired=#{paired ? "T" : "F"})"
      t = @r.pull "ttest$statistic"
      p = @r.pull "ttest$p.value"
      if (1-significance_level > p)
        t
      else
        0
      end
    end
    
    def pvalue(array1, array2)
      @r.assign "v1",array1
      @r.assign "v2",array2
      @r.eval "ttest = t.test(as.numeric(v1),as.numeric(v2))"
      @r.pull "ttest$p.value"
    end
        
    
    def ttest_single_value(array1, value2, significance_level=0.95)
      @r.assign "v1",array1
      @r.eval "ttest = t.test(as.numeric(v1),conf.level=#{significance_level})"
      min = @r.pull "ttest$conf.int[1]"
      max = @r.pull "ttest$conf.int[2]"
      if value2 <= min
        1
      elsif value2 >= max
        -1
      else
        0
      end
    end
    
    
    # example: 
    # files = ["/tmp/box.svg","/tmp/box.png"]
    # data = [ [ :method, [4,4,5,5,4,3,2] ], [ :method2, [1,2,3,4,5,4,6] ], [ :asdf, [9,1,8,0,7,1,6] ] ]
    # boxplot(files, data, "comparison1" )
    #
    def boxplot(files, data, title="", hline=nil, param="")
      LOGGER.debug("r-util> create boxplot "+data.inspect)
      raise "no hashes, to keep order" if data.is_a?(Hash)
      max = -1
      min = 100000
      max_median = -1
      min_median = 100000
      max_median_idx = -1
      min_median_idx = -1
      data.size.times do |i|
        values = data[i][1]
        max = [max,values.size].max
        min = [min,values.size].min
        med = values.to_scale.median
        #puts "#{data[i][0]} median: #{med}"
        #puts data[i][1].inspect
        max_median = [max_median,med].max
        max_median_idx = i if max_median==med
        min_median = [min_median,med].min
        min_median_idx = i if min_median==med
        data[i] = [data[i][0]+"(#{values.size})",data[i][1]]
      end
      if min != max
        times = max/min.to_f
        raise "box-plot values do not have equal size #{min} <-> #{max}" if times.floor != times.ceil
        data.size.times do |i|
          m = data[i][0]
          values = data[i][1]
          data[i] = [ m, values*times.to_i ] if values.size<max
        end
        min = 100000
        data.each do |m,values|
          max = [max,values.size].max
          min = [min,values.size].min
        end
      end
      assign_dataframe("boxdata",data.collect{|e| e[1]}.transpose,nil,data.collect{|e| e[0].to_s})
      #@r.eval "print('median')"
      #data.size.times.each do |i|
      #  @r.eval "print(median(boxdata[,#{i+1}]))"
      #  @r.eval "print(boxdata[,#{i+1}])"
      #end 
      param_str = (param!=nil and param.size>0) ? ",#{param}" : ""
      hlines = []
      hlines << [hline,"'gray60'"] if hline
      hlines << [max_median,2+max_median_idx] 
      hlines << [min_median,2+min_median_idx]
      plot_to_files(files, hlines) do |file|
        @r.eval "superboxplot(boxdata,main='#{title}',col=rep(2:#{data.size+1})#{param_str})"
      end
    end

    # embedds feature values of two datasets into 2D and plots it
    #        
    def feature_value_plot(files, dataset_uri1, dataset_uri2, dataset_name1, dataset_name2,
        features=nil, subjectid=nil, waiting_task=nil, direct_plot=false, title=nil)
        
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
      
      if (direct_plot)
        raise unless features.size==2
        @r.eval "df.2d <- df"
      else
        @r.eval "df.2d <- plot_pre_process(df, method='sammon')"
      end
      waiting_task.progress(75) if waiting_task
      
      title = "Sammon embedding of #{features.size} features" unless title
      LOGGER.debug("r-util> - plot data")
      plot_to_files(files) do |file|
        @r.eval "plot_split( df.2d, split, names, main='#{title}',xlab='x',ylab='y')"
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
    def stratified_split( dataset, metadata={}, missing_values="NA", pct=0.3, subjectid=nil, seed=42, split_features=nil, anti_stratification=false )
      stratified_split_internal( dataset, metadata, missing_values, nil, pct, subjectid, seed, split_features, anti_stratification )
    end
    
    # stratified splits a dataset into k datasets according the feature values
    # all features are taken into account unless <split_features> is given
    # returns two arrays of datasets
    def stratified_k_fold_split( dataset, metadata={}, missing_values="NA", num_folds=10, subjectid=nil, seed=42, split_features=nil )
      stratified_split_internal( dataset, metadata, missing_values, num_folds, nil, subjectid, seed, split_features )
    end    
    
    private
    def stratified_split_internal( dataset, metadata={}, missing_values="NA", num_folds=nil, pct=nil, subjectid=nil, seed=42, split_features=nil, stratification="super" )
      raise "internal error" if num_folds!=nil and pct!=nil
      k_fold_split = num_folds!=nil
      if k_fold_split
        raise "num_folds not a fixnum: #{num_folds}" unless num_folds.is_a?(Fixnum)
      else
        raise "pct is not a numeric: #{pct}" unless pct.is_a?(Numeric)
      end
      raise "not a loaded ot-dataset" unless dataset.is_a?(OpenTox::Dataset) and dataset.compounds.size>0 and dataset.features.size>0
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
        raise unless stratification=~/^(super|super4|anti)$/
        anti = ""
        super_method = ""
        preprocess = ""
        case stratification
        when "anti"
          anti = "anti_"
        when "super"
          super_method = ", method='cluster_knn'"
        when "super4"
          super_method = ", method='cluster_hierarchical'"
          preprocess = ", preprocess='pca'"
        end
        puts "split <- #{anti}stratified_split(#{df}, ratio=#{pct}, #{str_split_features} #{super_method} #{preprocess})"
        @r.eval "split <- #{anti}stratified_split(#{df}, ratio=#{pct}, #{str_split_features} #{super_method} #{preprocess})"
        split = @r.pull 'split'
        metadata[DC.title] = "Training dataset split of "+dataset.uri
        train = split_to_dataset( df, split, metadata, subjectid ){ |i| i==1 }
        metadata[DC.title] = "Test dataset split of "+dataset.uri
        test = split_to_dataset( df, split, metadata, subjectid, missing_values ){ |i| i==0 }
        return train, test
      end
    end
    public
    
    # dataset should be loaded completely (use Dataset.find)
    # takes duplicates into account
    # replaces missing values with param <missing_value>
    # returns dataframe-variable-name in R
    def dataset_to_dataframe( dataset, missing_values="NA", subjectid=nil, features=nil )
      LOGGER.debug "r-util> convert dataset to dataframe #{dataset.uri}"
      
      # count duplicates
      num_compounds = {}
      dataset.features.keys.each do |f|
        dataset.compounds.each do |c|
          if dataset.data_entries[c]
            val = dataset.data_entries[c][f]
            size = val==nil ? 1 : val.size
            num_compounds[c] = num_compounds[c]==nil ? size : [num_compounds[c],size].max
          else
            num_compounds[c] = 1
          end
        end
      end  
      
      # use either all, or the provided features, sorting is important as col-index := features
      if features
        features.sort!
      else
        features = dataset.features.keys.sort
      end
      compounds = []
      compound_names = []
      dataset.compounds.each do |c|
        count = 0
        num_compounds[c].times do |i|
          compounds << c
          compound_names << "#{c}$#{count}"
          count+=1
        end
      end

      #LOGGER.debug "converting to array"
      
      # values into 2D array, then to dataframe
      d_values = []
      dataset.compounds.each do |c|
        num_compounds[c].times do |i|
          c_values = []
          features.each do |f|
            if dataset.data_entries[c]
              val = dataset.data_entries[c][f]
              v = val==nil ? "" : val[i].to_s
            else
              raise "wtf" if i>0
              v = ""
            end
            v = missing_values if v.size()==0
            c_values << v
          end
          d_values << c_values
        end
      end  
      
      #LOGGER.debug "assigning"
      
      df_name = "df_#{dataset.uri.split("/")[-1].split("?")[0]}"
      assign_dataframe(df_name,d_values,compound_names,features)
      
      #LOGGER.debug "setting types"
      
      # set dataframe column types accordingly
      f_count = 1 #R starts at 1
      features.each do |f|
        if f=~/\/feature\/bbrc\//
          numeric=true
        else
          type = dataset.features[f][RDF.type]
          unless type
            LOGGER.debug "r-util> derive feature type by rest-call"
            feat = OpenTox::Feature.find(f,subjectid)
            type = feat.metadata[RDF.type]
          end
          numeric = type.to_a.flatten.include?(OT.NumericFeature)
        end
        unless numeric
          @r.eval "#{df_name}[,#{f_count}] <- as.character(#{df_name}[,#{f_count}])"
        else
          @r.eval "#{df_name}[,#{f_count}] <- as.numeric(#{df_name}[,#{f_count}])"
        end
        f_count += 1
      end
      #@r.eval "head(#{df_name})"
      #@r.eval "save.image(\"/tmp/image.R\")"
      
      # store compounds, and features (including metainformation)
      @@feats[df_name] = {}
      features.each do |f|
        @@feats[df_name][f] = dataset.features[f]
      end
      df_name
    end
    
    # converts a dataframe into a dataset (a new dataset is created at the dataset webservice)
    # this is only possible if a superset of the dataframe was created by dataset_to_dataframe (metadata and URIs!)
    def dataframe_to_dataset( df, metadata={}, subjectid=nil, missing_values="NA" )
      dataframe_to_dataset_indices( df, metadata, subjectid, nil, missing_values )
    end
    
    NEW = false
    
    private
    def dataframe_to_dataset_indices( df, metadata={}, subjectid=nil, compound_indices=nil, missing_values="NA" )
      raise unless @@feats[df].size>0

      missing_value_regexp = Regexp.new("^#{missing_values.to_s=="0" ? "(0.0|0)" : missing_values.to_s}$") unless NEW
      values, compound_names, features = pull_dataframe(df,missing_values)
      compounds = compound_names.collect{|c| c.split("$")[0]}
      features.each{|f| raise unless @@feats[df][f]}
      dataset = OpenTox::Dataset.create(CONFIG[:services]["opentox-dataset"],subjectid)
      dataset.add_metadata(metadata)
      LOGGER.debug "r-util> convert dataframe to dataset #{dataset.uri}"
      compounds.size.times{|i| dataset.add_compound(compounds[i]) if compound_indices==nil or compound_indices.include?(i)}
      features.each{|f| dataset.add_feature(f,@@feats[df][f])}
      features.size.times do |c|
        LOGGER.debug "r-util> dataframe to dataset - feature #{c+1} / #{features.size}" if 
          c%25==0 && (features.size*compounds.size)>100000
        if features[c]=~/\/feature\/bbrc\//
          numeric=true
        else
          type = @@feats[df][features[c]][RDF.type]
          unless type
            LOGGER.debug "r-util> derive feature type by rest-call"
            feat = OpenTox::Feature.find(features[c],subjectid)
            type = feat.metadata[RDF.type]
          end
          numeric = type.to_a.flatten.include?(OT.NumericFeature)
        end
        compounds.size.times do |r|
          if compound_indices==nil or compound_indices.include?(r)
            dataset.add(compounds[r],features[c],numeric ? values[r][c].to_f : values[r][c], true) if 
              ((NEW and values[r][c]!=nil) or (values[r][c]!="NA" and !(values[r][c] =~ missing_value_regexp)))
          end 
        end
      end
      dataset.save(subjectid)
      dataset
    end    
    
    def split_to_dataset( df, split, metadata={}, subjectid=nil, missing_values="NA" )
      indices = []
      split.size.times{|i| indices<<i if yield(split[i]) }
      dataset = dataframe_to_dataset_indices( df, metadata, subjectid, indices, missing_values )
      LOGGER.debug("r-util> split into #{dataset.uri}, c:#{dataset.compounds.size}, f:#{dataset.features.size}")
      dataset
    end
    
    def pull_dataframe(df,missing_values="NA")
      missing_value_regexp = Regexp.new("^#{missing_values.to_s=="0" ? "(0.0|0)" : missing_values.to_s}$") if NEW
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
           if NEW
             vals = row.chomp.gsub(missing_value_regexp,"").split("#").collect{|e| e.gsub("\"","")}
             compounds << vals[0]
             res << vals[1..-1].collect{|s| s=="" ? nil : s}
           else
             vals = row.chomp.split("#").collect{|e| e.gsub("\"","")}
             compounds << vals[0]
             res << vals[1..-1]
           end
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
    
    @@svg_plot_width = 12
    @@svg_plot_height = 8
        
    public
    def set_svg_plot_size(width,height)
      @@svg_plot_width = width
      @@svg_plot_height = height
    end
    
    private
    def plot_to_files(files,hlines=nil)
      files.each do |file|
        if file=~/(?i)\.svg/
          @r.eval("svg('#{file}',#{@@svg_plot_width},#{@@svg_plot_height})")
        elsif file=~/(?i)\.png/
          @r.eval("png('#{file}')")
        else
          raise "invalid format: "+file.to_s
        end
        yield file
        hlines.each{|hline,col| @r.eval("abline(h=#{hline}, col = #{col}, lty=2)")} if hlines
        LOGGER.debug "r-util> plotted to #{file}"
        @r.eval("dev.off()")
      end
    end
  end
end


