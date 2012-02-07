# pending: package dir hack ---------
# CONFIG[:base_dir] = "/home/<user>/opentox-ruby/www"
# PACKAGE_DIR = "/home/<user>/opentox-ruby/r-packages"
package_dir = CONFIG[:base_dir].split("/")
package_dir[-1] = "r-packages"
package_dir = package_dir.join("/")
PACKAGE_DIR = package_dir

require "tempfile"

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
    # fast_plot = true -> PCA, fast_plot = false -> SMACOF (iterative optimisation method) 
    #        
    def feature_value_plot(files, dataset_uri1, dataset_uri2, dataset_name1, dataset_name2,
        features=nil, fast_plot=true, subjectid=nil, waiting_task=nil)
        
      raise "r-package smacof missing" if fast_plot==false and !package_installed?("smacof")
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
      @r.eval "df.2d <- plot_pre_process(df, method='#{(fast_plot ? "pca" : "smacof")}')"
      waiting_task.progress(75) if waiting_task
      
      if fast_plot
        info = "main='PCA-Embedding of #{features.size} features',xlab='PC1',ylab='PC2'"
      else
        info = "main='SMACOF-Embedding of #{features.size} features',xlab='x',ylab='y'"
      end
      LOGGER.debug("r-util> - plot data")
      plot_to_files(files) do |file|
        @r.eval "plot_split( df.2d, split, names, #{info})"
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
    
    # stratified splits a dataset into two dataset the feature values
    # all features are taken into account unless <split_features> is given
    def stratified_split( dataset, missing_values="NA", pct=0.3, subjectid=nil, seed=42, split_features=nil )
      raise "not a loaded ot-dataset" unless dataset.is_a?(OpenTox::Dataset) and dataset.compounds.size>0 and dataset.features.size>0
      LOGGER.debug("r-util> apply stratified split to #{dataset.uri}")
      
      df = dataset_to_dataframe( dataset, missing_values, subjectid, split_features )
      @r.eval "set.seed(#{seed})"
      @r.eval "split <- stratified_split(#{df}, ratio=#{pct})"
      split = @r.pull 'split'
      split = split.collect{|s| 1-s.to_i} # reverse 1s and 0s, as 1 means selected, but 0 will be first set
      split_to_datasets( df, split, subjectid )
    end
    
    # dataset should be loaded completely (use Dataset.find)
    # takes duplicates into account
    # replaces missing values with param <missing_value>
    # returns dataframe-variable-name in R
    def dataset_to_dataframe( dataset, missing_value="NA", subjectid=nil, features=nil )
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
      dataset.compounds.each do |c|
        num_compounds[c].times do |i|
          compounds << c
        end
      end

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
            v = missing_value if v.size()==0
            c_values << v
          end
          d_values << c_values
        end
      end  
      df_name = "df_#{dataset.uri.split("/")[-1].split("?")[0]}"
      assign_dataframe(df_name,d_values,compounds,features)
      
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
    def dataframe_to_dataset( df, subjectid=nil )
      dataframe_to_dataset_indices( df, subjectid, nil)
    end
    
    private
    def dataframe_to_dataset_indices( df, subjectid=nil, compound_indices=nil )
      raise unless @@feats[df].size>0
      values, compounds, features = pull_dataframe(df)
      features.each{|f| raise unless @@feats[df][f]}
      dataset = OpenTox::Dataset.create(CONFIG[:services]["opentox-dataset"],subjectid)
      LOGGER.debug "r-util> convert dataframe to dataset #{dataset.uri}"
      compounds.size.times{|i| dataset.add_compound(compounds[i]) if compound_indices==nil or compound_indices.include?(i)}
      features.each{|f| dataset.add_feature(f,@@feats[df][f])}
      features.size.times do |c|
        feat = OpenTox::Feature.find(features[c],subjectid)
        nominal = feat.metadata[RDF.type].to_a.flatten.include?(OT.NominalFeature)
        compounds.size.times do |r|
          if compound_indices==nil or compound_indices.include?(r)
            dataset.add(compounds[r],features[c],nominal ? values[r][c] : values[r][c].to_f) if values[r][c]!="NA"
          end 
        end
      end
      dataset.save(subjectid)
      dataset
    end    
    
    def split_to_datasets( df, split, subjectid=nil )
      sets = []
      (split.min.to_i .. split.max.to_i).each do |i|
        indices = []
        split.size.times{|j| indices<<j if split[j]==i}
        dataset = dataframe_to_dataset_indices( df, subjectid, indices )
        LOGGER.debug("r-util> split into #{dataset.uri}, c:#{dataset.compounds.size}, f:#{dataset.features.size}")
        sets << dataset
      end
      sets
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


