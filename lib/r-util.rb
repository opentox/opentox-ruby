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

class RinRuby
  def puts(object)
    object.to_s.split("\n").each do |s|
      LOGGER.debug "R> #{s.chomp}" if s.chomp.length>0
    end
  end
end

module OpenTox
  
  class RUtil
    
    @@feats = {}
      
    def initialize
      @r = RinRuby.new(true,false) unless defined?(@r) and @r
      #@r.eval "sink(type='message')"
      @r.eval ".libPaths('#{PACKAGE_DIR}')"
      @r_packages = @r.pull "installed.packages()[,1]"
      ["sampling","gam","vegan","dynamicTreeCut","proxy"].each{|l| install_package(l)} #"caret", "smacof", "TunePareto"
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
    
#    def ttest_matrix_deviation(arrays, significance_level=0.95)
#    
#      LOGGER.info("perform ttest matrix deviation")
#      result = []
#      arrays.size.times do |i|
#        result[i] = [] 
#        @r.assign "v#{i}",arrays[i]
#        @r.eval "v#{i}<-abs(as.numeric(v#{i}))"
#      end 
#      (arrays.size-1).times do |i|
#        (i+1..arrays.size-1).each do |j|
#          raise if arrays[i].size!=arrays[j].size
#          @r.eval "ttest = t.test(v#{i}-v#{j},conf.level=#{significance_level})"
#          min = @r.pull "ttest$conf.int[1]"
#          max = @r.pull "ttest$conf.int[2]"
#          if (min>0 or max<0)
#            result[i][j]=true
#            result[j][i]=true
#          else
#            result[i][j]=false
#            result[j][i]=false
#          end            
#        end
#      end
#      result
#    end
    
    def pvalue_test_matrix(test, arrays, significance_level=0.95, params="")
      LOGGER.info("perform test '#{test}' matrix")
      result = []
      arrays.size.times do |i|
        result[i] = [] 
        @r.assign "v#{i}",arrays[i]
        @r.eval "v#{i}<-as.numeric(v#{i})"
      end 
      (arrays.size-1).times do |i|
        (i+1..arrays.size-1).each do |j|
          @r.eval "test = #{test}(v#{i},v#{j},#{params})"
          t = @r.pull "test$statistic"
          p = @r.pull "test$p.value"
          if (1-significance_level > p)
            result[i][j]=true
            result[j][i]=true
          else
            result[i][j]=false
            result[j][i]=false
          end            
        end
      end
      result
    end
    
    def ttest_matrix(arrays, paired, significance_level=0.95)
      (arrays.size-1).times do |i|
        (i+1..arrays.size-1).each do |j|
          raise if paired && arrays[i].size!=arrays[j].size
        end   
      end           
      params = "paired=#{paired ? "T" : "F"}"
      pvalue_test_matrix("t.test",arrays,significance_level,params)
    end    
    
    def ftest_matrix(arrays, paired, significance_level=0.95)
      pvalue_test_matrix("var.test",arrays,significance_level)
    end  
    
    def ttest_closer_to_zero_matrix(arrays, paired, significance_level=0.95)
      (arrays.size-1).times do |i|
        (i+1..arrays.size-1).each do |j|
          raise if paired && arrays[i].size!=arrays[j].size
        end   
      end           
      params = "paired=#{paired ? "T" : "F"}"
      pvalue_test_matrix("ttest_closer_to_zero",arrays,significance_level,params)
    end      
    
    def pvalue_test(test, array1, array2, significance_level=0.95, params="")
      LOGGER.info("perform test '#{test}'")
        @r.assign "v1",array1
        @r.assign "v2",array2
        @r.eval "test = #{test}(as.numeric(v1),as.numeric(v2),#{params})"
        t = @r.pull "test$statistic"
        p = @r.pull "test$p.value"
        if (1-significance_level > p)
          t
        else
          0
        end
    end
    
    # <0 -> array1 << array2
    # 0  -> no significant difference
    # >0 -> array2 >> array1
    def ttest(array1, array2, paired, significance_level=0.95)
      raise if paired && array1.size!=array2.size
      params = "paired=#{paired ? "T" : "F"}"
      pvalue_test("t.test",array1,array2,significance_level,params)
    end
    
    def ftest(array1, array2, significance_level=0.95)
      pvalue_test("var.test",array1,array2,significance_level)
    end
    
    def ttest_closer_to_zero
      raise if paired && array1.size!=array2.size
      params = "paired=#{paired ? "T" : "F"}"
      pvalue_test("ttest_closer_to_zero",array1,array2,significance_level,params)
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
        LOGGER.debug "perform ttest-single: significant=true, #{value2} is lower than conf-interval [ #{min} - #{max} ]"
        1
      elsif value2 >= max
        LOGGER.debug "perform ttest-single: significant=true, #{value2} is higher than conf-interval [ #{min} - #{max} ]"
        -1
      else
        LOGGER.debug "perform ttest-single: significant=false, #{value2} is inside conf-interval [ #{min} - #{max} ]"
        0
      end
    end
    
    
    private
    def get_r_cols(pair_colors=false)
      cols = ["red","cyan","green","magenta","blue","orange","seagreen","salmon","goldenrod","gray","orchid","khaki"]
      if pair_colors
        pair_cols=[]
        cols.each{|c| pair_cols<<c; pair_cols<<"dark#{c}"}
        cols = pair_cols
      end
      "col=c('#{cols.join("','")}')"
    end
    
    public
    # example: 
    # files = ["/tmp/box.svg","/tmp/box.png"]
    # data = [ [ :method, [4,4,5,5,4,3,2] ], [ :method2, [1,2,3,4,5,4,6] ], [ :asdf, [9,1,8,0,7,1,6] ] ]
    # boxplot(files, data, "comparison1" )
    #
    def boxplot(files, data, title="", hline=nil, param="", pair_colors=false)
      LOGGER.debug("r-util> create boxplot "+data.inspect)
      raise "no hashes, to keep order" if data.is_a?(Hash)
      raise "boxplot: data is empty" if data.size==0
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
        data[i] = [data[i][0].to_s+"(#{values.size})",data[i][1]] if @@boxplot_alg_info
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
        #@r.eval "superboxplot(boxdata,alg_info=#{@@boxplot_alg_info ? "T" : "F"},main='#{title}',col=rep(2:#{data.size+1})#{param_str})"
        @r.eval "superboxplot(boxdata,alg_info=#{@@boxplot_alg_info ? "T" : "F"},main='#{title}',#{get_r_cols(pair_colors)}#{param_str})"
      end
    end

    # embedds feature values of two datasets into 2D and plots it
    #        
    def feature_value_plot(files, dataset_uri1, dataset_uri2, dataset_name1, dataset_name2, feature_type,
        prediction_feature=nil, subjectid=nil, waiting_task=nil, direct_plot=false, title=nil, color_feature=nil )
        
      raise "feature_type has to be binary or numerical" unless ["binary","numerical"].include?(feature_type)
        
      LOGGER.debug("r-util> create feature value plot #{feature_type}")
      d1 = OpenTox::Dataset.find(dataset_uri1,subjectid)
      d2 = OpenTox::Dataset.find(dataset_uri2,subjectid)
      
      raise "different\n#{d1.features.keys.sort.to_yaml}\n#{d2.features.keys.sort.to_yaml}" if 
          (d1.features.keys.sort != d2.features.keys.sort)
      features = d1.features.keys
      if prediction_feature
        if features.include?(prediction_feature)
          features -= [prediction_feature]
        else
          LOGGER.debug "prediction feature #{prediction_feature} cannot be remvoed because not included in #{dataset_uri1}"
        end
      end 
      
      raise "at least two features needed" if d1.features.keys.size<2
      waiting_task.progress(25) if waiting_task
      
      df1 = dataset_to_dataframe(d1,0,subjectid,features)
      df2 = dataset_to_dataframe(d2,0,subjectid,features)
      waiting_task.progress(50) if waiting_task
      
      @r.eval "df <- rbind(#{df1},#{df2})"
      @r.eval "split <- c(rep(1,nrow(#{df1})),rep(0,nrow(#{df2})))"
      @r.names = [dataset_name1, dataset_name2]
      LOGGER.debug("r-util> - convert data to 2d")
      #@r.eval "save.image(\"/tmp/image.R\")"
      
      if (color_feature)
        color = []
        [d1,d2].each do |d|
          raise "no #{color_feature}, instead: #{d.features.keys.sort.inspect}" unless d.features.has_key?(color_feature)
          d.compounds.each do |c|
            color += d.data_entries[c][color_feature]
          end
        end
        @r.assign "color",color
      end
      
      if (direct_plot)
        raise unless features.size==2
        @r.eval "df.2d <- df"
        x = features[0].split("/")[-1]
        y = features[1].split("/")[-1]
      else
        @r.eval "df.2d <- plot_pre_process(df, '#{feature_type}', method='sammon')"
        x = "x"
        y = "y"
      end
      waiting_task.progress(75) if waiting_task
      
      title = "Sammon embedding of #{features.size} features" unless title
      LOGGER.debug("r-util> - plot data")
      plot_to_files(files) do |file|
        if (color_feature)
          @r.eval "plot_split( df.2d, color_idx=color, circle_idx=split, main='#{title}',xlab='#{x}',ylab='#{y}')"
        else
          @r.eval "plot_split( df.2d, color_idx=split, main='#{title}',xlab='#{x}',ylab='#{y}')"
        end
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
          if( log && ( min(data1)<=0 || min(data2)<=0 ))
          {
            print('disabling log because of datapoints <= 0')
            log = FALSE
          }   
          if (log)
          {
            data1 <- log(data1)
            data2 <- log(data2)
            xlab = paste('logarithm of',xlab,sep=' ')
          }
          xlims <- round(c(min(c(min(data1),min(data2))),max(c(max(data1),max(data2)))))
          save.image('/tmp/image.R')
          h <- hist(c(data1,data2),plot=F)
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
    def stratified_split( dataset, metadata={}, missing_values="NA", pct=0.3, subjectid=nil, 
      seed=42, split_features=nil, anti_stratification=false, store_split_clusters=false )
      stratified_split_internal( dataset, metadata, missing_values, nil, pct, subjectid, 
        seed, split_features, anti_stratification, store_split_clusters )
    end
    
    # stratified splits a dataset into k datasets according the feature values
    # all features are taken into account unless <split_features> is given
    # returns two arrays of datasets
    def stratified_k_fold_split( dataset, metadata={}, missing_values="NA", num_folds=10, subjectid=nil, seed=42, split_features=nil )
      stratified_split_internal( dataset, metadata, missing_values, num_folds, nil, subjectid, seed, split_features )
    end    
    
    private
    def stratified_split_internal( dataset, metadata={}, missing_values="NA", num_folds=nil, 
        pct=nil, subjectid=nil, seed=42, split_features=nil, stratification="super", store_split_clusters=false )
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
        raise unless stratification=~/^(super|super4|super5|super_bin|contra_eucl2|contra_bin2)$/
        anti = ""
        super_method = ""
        super_method_2 = ""
        #preprocess = ""
        case stratification
        when "contra_eucl2"
          feature_type = "numerical"
          anti = "contra_"
        when "contra_bin2"
          feature_type = "binary"
          anti = "contra_"
        when "super"
          feature_type = "numerical"
          super_method = ", method='cluster_knn'"
        when "super4"
          feature_type = "numerical"
          super_method = ", method='cluster_hierarchical'"
        when "super5"
          feature_type = "numerical"
          super_method = ", method='cluster_hierarchical'"
          super_method_2 = ", method_2='explicit'"
        when "super_bin"
          feature_type = "binary"
          super_method = ", method='cluster_hierarchical'"
          super_method_2 = ", method_2='explicit'"
        else
          raise "strat unknown"
        end
        cmd = "split <- #{anti}stratified_split(#{df}, '#{feature_type}', ratio=#{pct}, #{str_split_features} #{super_method} #{super_method_2})" # #{preprocess}
        LOGGER.debug cmd
        @r.eval cmd
        split = @r.pull 'split$split'
        cluster = (store_split_clusters ?  @r.pull('split$cluster') : nil)
        metadata[DC.title] = "Training dataset split of "+dataset.uri
        train = split_to_dataset( df, split, metadata, subjectid, missing_values, cluster ){ |i| i==1 }
        metadata[DC.title] = "Test dataset split of "+dataset.uri
        test = split_to_dataset( df, split, metadata, subjectid, missing_values, cluster ){ |i| i==0 }

        #f = "/tmp/split_pic.svg"
        #LOGGER.debug "plotting to #{f} .."
        #@r.eval "num_feats = #{split_features ? split_features.size : "ncol(plot_data)"}"
        #@r.eval "plot_data = process_data(#{df}, #{str_split_features})"
        #@r.eval "plot_data = plot_pre_process(plot_data, method='sammon')"
        #@r.eval "title = paste('sammon embedding for splitting #{df},',num_feats,'features,',nrow(plot_data),'instances')"
        #plot_to_files([f]) do |file|
        #  @r.eval "plot_split(plot_data,color_idx=split$split, main=title)"
        #end
        #LOGGER.debug "plotting to #{f} .. done"
          
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
    def dataframe_to_dataset_indices( df, metadata={}, subjectid=nil, compound_indices=nil, missing_values="NA", cluster=nil )
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
      features.size.times do |f_i|
        LOGGER.debug "r-util> dataframe to dataset - feature #{f_i+1} / #{features.size}" if 
          f_i%25==0 && (features.size*compounds.size)>100000
        if features[f_i]=~/\/feature\/bbrc\//
          numeric="int"
        else
          type = @@feats[df][features[f_i]][RDF.type]
          unless type
            LOGGER.debug "r-util> derive feature type by rest-call"
            feat = OpenTox::Feature.find(features[f_i],subjectid)
            type = feat.metadata[RDF.type]
          end
          numeric = type.to_a.flatten.include?(OT.NumericFeature) ? "float" : nil
        end
        case numeric
        when "int"
          def convert_numeric(v); v.to_i; end
        when "float"             
          def convert_numeric(v); v.to_f; end
        else
          def convert_numeric(v); v; end
        end
        compounds.size.times do |c_i|
          if compound_indices==nil or compound_indices.include?(c_i)
            dataset.add(compounds[c_i],features[f_i],convert_numeric(values[c_i][f_i]), true) if 
              ((NEW and values[c_i][f_i]!=nil) or (values[c_i][f_i]!="NA" and !(values[c_i][f_i] =~ missing_value_regexp)))
          end 
        end
        
      end
      if cluster
        cluster_feature = "http://no.such.domain/feature/split_cluster"
        dataset.add_feature(cluster_feature)
        #LOGGER.warn "adding feature #{cluster_feature}" 
        compounds.size.times do |r|
          if compound_indices==nil or compound_indices.include?(r)
            dataset.add(compounds[r],cluster_feature,cluster[r],true)
          end 
        end
      else
        #LOGGER.warn "no cluster feature" 
      end
      dataset.save(subjectid)
      dataset
    end    
    
    def split_to_dataset( df, split, metadata={}, subjectid=nil, missing_values="NA", cluster=nil )
      indices = []
      split.size.times{|i| indices<<i if yield(split[i]) }
      dataset = dataframe_to_dataset_indices( df, metadata, subjectid, indices, missing_values, cluster )
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
    
    @@png_plot_width = 800
    @@png_plot_height = 600
    @@png_plot_pointsize = 12
    
    def set_png_plot_size(width,height,pointsize)
      @@png_plot_width = width
      @@png_plot_height = height
      @@png_plot_pointsize = pointsize
    end
    
    @@boxplot_alg_info = true
    
    def set_boxplot_alg_info(boxplot_alg_info)
      @@boxplot_alg_info = boxplot_alg_info
    end
    
    private
    def plot_to_files(files,hlines=nil)
      files.each do |file|
        if file=~/(?i)\.svg/
          @r.eval("svg('#{file}',#{@@svg_plot_width},#{@@svg_plot_height})")
        elsif file=~/(?i)\.png/
          @r.eval("png('#{file}',width=#{@@png_plot_width},height=#{@@png_plot_height},pointsize=#{@@png_plot_pointsize})")
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


