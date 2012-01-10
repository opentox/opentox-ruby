# R integration
# workaround to initialize R non-interactively (former rinruby versions did this by default)
# avoids compiling R with X
R = nil
require "rinruby" 
require "statsample"
require 'uri'
require 'transform.rb'
require 'utils.rb'

module OpenTox

  # Wrapper for OpenTox Algorithms
  module Algorithm 

    include OpenTox

    # Execute algorithm with parameters, please consult the OpenTox API and the webservice documentation for acceptable parameters
    # @param [optional,Hash] params Algorithm parameters
    # @param [optional,OpenTox::Task] waiting_task (can be a OpenTox::Subtask as well), progress is updated accordingly
    # @return [String] URI of new resource (dataset, model, ...)
    def run(params=nil, waiting_task=nil)
      LOGGER.info "Running algorithm '"+@uri.to_s+"' with params: "+params.inspect
      RestClientWrapper.post(@uri, params, {:accept => 'text/uri-list'}, waiting_task).to_s
    end
    
    # Get OWL-DL representation in RDF/XML format
    # @return [application/rdf+xml] RDF/XML representation
    def to_rdfxml
      s = Serializer::Owl.new
      s.add_algorithm(@uri,@metadata)
      s.to_rdfxml
    end

    # Generic Algorithm class, should work with all OpenTox webservices
    class Generic 
      include Algorithm
      
      # Find Generic Opentox Algorithm via URI, and loads metadata, could raise NotFound/NotAuthorized error
      # @param [String] uri Algorithm URI
      # @return [OpenTox::Algorithm::Generic] Algorithm instance
      def self.find(uri, subjectid=nil)
        return nil unless uri
        alg = Generic.new(uri)
        alg.load_metadata( subjectid )
        raise "cannot load algorithm metadata" if alg.metadata==nil or alg.metadata.size==0
        alg
      end
      
    end

    # Fminer algorithms (https://github.com/amaunz/fminer2)
    class Fminer
      include Algorithm
      attr_accessor :prediction_feature, :training_dataset, :minfreq, :compounds, :db_class_sizes, :all_activities, :smi
      
      def check_params(params,per_mil,subjectid=nil)
        raise OpenTox::NotFoundError.new "Please submit a dataset_uri." unless params[:dataset_uri] and  !params[:dataset_uri].nil?
        raise OpenTox::NotFoundError.new "Please submit a prediction_feature." unless params[:prediction_feature] and  !params[:prediction_feature].nil?
        @prediction_feature = OpenTox::Feature.find params[:prediction_feature], subjectid
        @training_dataset = OpenTox::Dataset.find "#{params[:dataset_uri]}", subjectid
        raise OpenTox::NotFoundError.new "No feature #{params[:prediction_feature]} in dataset #{params[:dataset_uri]}" unless @training_dataset.features and @training_dataset.features.include?(params[:prediction_feature])

        unless params[:min_frequency].nil? 
          @minfreq=params[:min_frequency].to_i
          raise "Minimum frequency must be a number >0!" unless @minfreq>0
        else
          @minfreq=OpenTox::Algorithm.min_frequency(@training_dataset,per_mil) # AM sugg. 8-10 per mil for BBRC, 50 per mil for LAST
        end
      end

      def add_fminer_data(fminer_instance, params, value_map)

        id = 1 # fminer start id is not 0
        @training_dataset.data_entries.each do |compound,entry|
          begin
            smiles = OpenTox::Compound.smiles(compound.to_s)
          rescue
            LOGGER.warn "No resource for #{compound.to_s}"
            next
          end
          if smiles == '' or smiles.nil?
            LOGGER.warn "Cannot find smiles for #{compound.to_s}."
            next
          end
          
          value_map=params[:value_map] unless params[:value_map].nil?
          entry.each do |feature,values|
            if feature == @prediction_feature.uri
              values.each do |value|
                if value.nil? 
                  LOGGER.warn "No #{feature} activity for #{compound.to_s}."
                else
                  if @prediction_feature.feature_type == "classification"
                    activity= value_map.invert[value].to_i # activities are mapped to 1..n
                    @db_class_sizes[activity-1].nil? ? @db_class_sizes[activity-1]=1 : @db_class_sizes[activity-1]+=1 # AM effect
                  elsif @prediction_feature.feature_type == "regression"
                    activity= value.to_f 
                  end
                  begin
                    fminer_instance.AddCompound(smiles,id)
                    fminer_instance.AddActivity(activity, id)
                    @all_activities[id]=activity # DV: insert global information
                    @compounds[id] = compound
                    @smi[id] = smiles
                    id += 1
                  rescue Exception => e
                    LOGGER.warn "Could not add " + smiles + "\t" + value.to_s + " to fminer"
                    LOGGER.warn e.backtrace
                  end
                end
              end
            end
          end
        end
      end

    end

      # Backbone Refinement Class mining (http://bbrc.maunz.de/)
      class BBRC < Fminer
        # Initialize bbrc algorithm
        def initialize(subjectid=nil)
          super File.join(CONFIG[:services]["opentox-algorithm"], "fminer/bbrc")
          load_metadata(subjectid)
        end
      end

      # LAtent STructure Pattern Mining (http://last-pm.maunz.de)
      class LAST < Fminer
        # Initialize last algorithm
        def initialize(subjectid=nil)
          super File.join(CONFIG[:services]["opentox-algorithm"], "fminer/last")
          load_metadata(subjectid)
        end
      end


    # Create lazar prediction model
    class Lazar
      include Algorithm
      # Initialize lazar algorithm
      def initialize(subjectid=nil)
        super File.join(CONFIG[:services]["opentox-algorithm"], "lazar")
        load_metadata(subjectid)
      end
    end


    # Structural Graph Clustering by TU Munich
    # Finds clusters similar to a query structure in a given training dataset
    # May be queried for cluster membership of an unknown compound
    class StructuralClustering
      attr_accessor :training_dataset_uri, :training_threshold, :query_dataset_uri, :query_threshold, :target_clusters_array

      # @params[String] Training dataset_uri
      # @params[Float]  Similarity threshold for training (optional)
      # @params[String] Cluster service uri (no AA)
      def initialize training_dataset_uri, training_threshold=0.8, cluster_service_uri = "http://opentox-dev.informatik.tu-muenchen.de:8080/OpenTox/algorithm/StructuralClustering"

        if (training_dataset_uri =~ URI::regexp).nil? || (cluster_service_uri =~ URI::regexp).nil? 
          raise "Invalid URI."
        end
        @training_dataset_uri = training_dataset_uri
        if !OpenTox::Algorithm.numeric? training_threshold || training_threshold <0 || training_threshold >1
          raise "Training threshold out of bounds."
        end
        @training_threshold = training_threshold.to_f

        # Train a cluster model
        params = {:dataset_uri => @training_dataset_uri, :threshold => @training_threshold }
        @cluster_model_uri = OpenTox::RestClientWrapper.post cluster_service_uri, params
        cluster_model_rdf = OpenTox::RestClientWrapper.get @cluster_model_uri
        @datasets = OpenTox::Parser::Owl.from_rdf cluster_model_rdf, OT.Dataset, true # must extract OT.Datasets from model

        # Process parsed OWL objects
        @clusterid_dataset_map = Hash.new
        @datasets.each { |d|
          begin
            d.metadata[OT.hasSource]["Structural Clustering cluster "] = "" # must parse in metadata for string (not elegant)
            @clusterid_dataset_map[d.metadata[OT.hasSource].to_i] = d.uri
          rescue Exception => e
            # ignore other entries!
          end
        }
      end

      # Whether a model has been trained
      def trained?
        !@cluster_model_uri.nil?
      end

      # Instance query: clusters for a compound
      # @params[String] Query compound
      # @params[Float]  Similarity threshold for query to clusters (optional)
      def get_clusters query_compound_uri, query_threshold = 0.5

        if !OpenTox::Algorithm.numeric? query_threshold || query_threshold <0 || query_threshold >1
          raise "Query threshold out of bounds."
        end
        @query_threshold = query_threshold.to_f


        # Preparing a query dataset
        query_dataset = OpenTox::Dataset.new
        @query_dataset_uri = query_dataset.save
        query_dataset = OpenTox::Dataset.find @query_dataset_uri
        query_dataset.add_compound query_compound_uri
        @query_dataset_uri = query_dataset.save

        # Obtaining a clustering for query compound
        params = { :dataset_uri => @query_dataset_uri, :threshold => @query_threshold }
        cluster_query_dataset_uri = OpenTox::RestClientWrapper.post @cluster_model_uri, params
        cluster_query_dataset = OpenTox::Dataset.new cluster_query_dataset_uri
        cluster_query_dataset.load_all

        # Reading cluster ids for features from metadata
        feature_clusterid_map = Hash.new
        pattern="Prediction feature for cluster assignment " # must parse for string in metadata (not elegant)
        cluster_query_dataset.features.each { |feature_uri,metadata|
          metadata[DC.title][pattern]=""
          feature_clusterid_map[feature_uri] = metadata[DC.title].to_i
        }
        
        # Integrity check
        unless cluster_query_dataset.compounds.size == 1
          raise "Number of predicted compounds is != 1."
        end

        # Process data entry
        query_compound_uri = cluster_query_dataset.compounds[0]
        @target_clusters_array = Array.new
        cluster_query_dataset.features.keys.each { |cluster_membership_feature|
        
          # Getting dataset URI for cluster
          target_cluster = feature_clusterid_map[cluster_membership_feature]
          dataset = @clusterid_dataset_map[target_cluster]
        
          # Finally look up presence
          data_entry = cluster_query_dataset.data_entries[query_compound_uri]
          present = data_entry[cluster_membership_feature][0]

          # Store result
          @target_clusters_array << dataset if present > 0.5 # 0.0 for absence, 1.0 for presence
        }
      end

    end



    module Neighbors

      # Local multi-linear regression (MLR) prediction from neighbors.
      # Uses propositionalized setting.
      # @param [Hash] params Keys `:acts, :sims, :props` are required, `:maxcols` optional.
      # @return [Numeric] A prediction value.
      def self.local_mlr_prop(params)

        confidence=0.0
        prediction=nil

        LOGGER.debug "Local MLR."
        if params[:acts].size>0

          acts = params[:acts].collect
          sims = params[:sims][1].collect
          n_prop = params[:props][0].collect
          q_prop = params[:props][1].collect
          props = [ n_prop, q_prop ]
          maxcols = ( params[:maxcols].nil? ? (sims.size/3.0).ceil : params[:maxcols] )

          # Transform y
          acts_autoscaler = OpenTox::Transform::LogAutoScale.new(acts.to_gv)
          acts = acts_autoscaler.vs.to_a

          # Predict
          prediction = pcr( {:n_prop => props[0], :q_prop => props[1], :sims => sims, :acts => acts, :maxcols => maxcols} )

          # Restore
          prediction = acts_autoscaler.restore( [ prediction ].to_gv )[0]
          prediction = nil if (!prediction.nil? && prediction.infinite?)
          LOGGER.debug "Prediction is: '" + prediction.to_s + "'."
          params[:conf_stdev] = false if params[:conf_stdev].nil?
          confidence = get_confidence({:sims => sims, :acts => acts, :conf_stdev => params[:conf_stdev]})
          confidence = 0.0 if prediction.nil?
        end

        {:prediction => prediction, :confidence => confidence}
      end

      # Multi-linear regression (unweighted).
      # Robust, multivariate (X) and univariate (y) outlier detection.
      # @param [Hash] params Keys `:n_prop, :q_prop, :sims, :acts, :maxcols` are required.
      # @return [Numeric] A prediction value.
      def self.mlr(params)

        # Uses Statsample Library (http://ruby-statsample.rubyforge.org/) by C. Bustos
        # Statsample operations build on GSL and offer an R-like access to data
        LOGGER.debug "MLR..."
        begin
          n_prop = params[:n_prop].collect
          q_prop = params[:q_prop].collect
          acts = params[:acts].collect
          sims = params[:sims].collect
          maxcols = params[:maxcols]

          nr_cases, nr_features = n_prop.size, n_prop[0].size
          maxcols = nr_features if maxcols > nr_features

          data_matrix = GSL::Matrix.alloc(n_prop.flatten, nr_cases, nr_features)
          query_matrix = GSL::Matrix.alloc(q_prop.flatten, 1, nr_features) # same nr_features

          @r = RinRuby.new(false,false)   # global R instance leads to Socket errors after a large number of requests
          outliers = OpenTox::Algorithm::Similarity.outliers( { :query_matrix => query_matrix, :data_matrix => data_matrix, :acts => acts, :r => @r  } )
          LOGGER.debug "Detected #{outliers.size} outliers: [#{outliers.join(", ")}]"
          #if (outliers.include?(-1))
          #  raise "Query is an outlier."
          #end

          temp_dm = []; temp_acts = []; temp_sims = []
          data_matrix.to_a.each_with_index { |elem, idx| temp_dm << elem unless outliers.include? idx }
          nr_cases, nr_features = temp_dm.size, temp_dm[0].size
          data_matrix = GSL::Matrix.alloc(temp_dm.flatten, nr_cases, nr_features)
          acts.each_with_index { |elem, idx| temp_acts << elem unless outliers.include? idx }
          acts = temp_acts # same nr_features
          sims.each_with_index { |elem, idx| temp_sims << elem unless outliers.include? idx }
          sims = temp_sims # same nr_features

          @r.eval 'fstr <- "y ~ ."'
          @r.x = data_matrix.to_a.flatten
          @r.y = acts.to_a.flatten
          @r.w = sims.to_a.flatten
          @r.q = query_matrix.to_a.flatten

          @r.eval "x <- matrix(x, #{nr_cases}, #{nr_features}, byrow=T)"
          @r.eval 'df <- data.frame(y,x)'

          # optimize selection features; adjust query accordingly
          begin
            LOGGER.debug "Best subset..."
            @r.eval 'suppressPackageStartupMessages(library("leaps"))'
            @r.eval "allss = summary( regsubsets( as.formula(fstr), data=df, nvmax=#{[ (nr_cases / 3).floor, nr_features ].min}, method=\"exhaustive\", weights=w) )"
            @r.eval 'idx = as.vector(allss$which[which.max(allss$adjr2),])'
            @r.eval 'idx[1] = T' # enforce intercept
            @r.eval 'intidx = as.integer(idx)'
          rescue Exception => e
            LOGGER.debug "#{e.class}: #{e.message}"
            LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
          end
          LOGGER.debug "Indices: [" + @r.intidx.to_a.join(", ") + "]"
          raise "No features left" if (@r.intidx.to_a.inject{|sum,elem| sum + elem}) == 1
          @r.eval 'suppressPackageStartupMessages(library("MASS"))'
          @r.eval 'df <- df[,idx]'
          @r.eval 'fit <- rlm( as.formula(fstr), data=df, psi = psi.bisquare, weights=w, wt.method="case")'
          @r.eval 'q <- q[idx[2:length(idx)]]'
          @r.eval 'q <- data.frame( matrix( q, 1, length(q) ) )'
          @r.eval 'names(q) = names(df)[2:length(names(df))]'
          @r.eval 'pred <- predict(fit, q, interval="confidence")'
          
          (@r.pred.to_a.flatten)[0] # [1] is lwr, [2] upr confidence limit.

        rescue Exception => e
          LOGGER.debug "#{e.class}: #{e.message}"
          #LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
        end

      end


      # Partial least squares regression.
      # Robust, multivariate (X) and univariate (y) outlier detection.
      # @param [Hash] params Keys `:n_prop, :q_prop, :sims, :acts, :maxcols` are required.
      # @return [Numeric] A prediction value.
      def self.pcr(params)

        # Uses Statsample Library (http://ruby-statsample.rubyforge.org/) by C. Bustos
        # Statsample operations build on GSL and offer an R-like access to data
        LOGGER.debug "PCR..."
        begin
          n_prop = params[:n_prop].collect
          q_prop = params[:q_prop].collect
          acts = params[:acts].collect
          maxcols = params[:maxcols]

          nr_cases, nr_features = n_prop.size, n_prop[0].size
          maxcols = nr_features if maxcols > nr_features

          data_matrix = GSL::Matrix.alloc(n_prop.flatten, nr_cases, nr_features)
          query_matrix = GSL::Matrix.alloc(q_prop.flatten, 1, nr_features) # same nr_features

          ### Model
          @r = RinRuby.new(false,false)   # global R instance leads to Socket errors after a large number of requests
          @r.eval 'suppressPackageStartupMessages(library("pls"))'
          outliers = OpenTox::Algorithm::Similarity.outliers( { :query_matrix => query_matrix, :data_matrix => data_matrix, :acts => acts, :r => @r  } )


          LOGGER.debug "Detected #{outliers.size} outliers: [#{outliers.join(", ")}]"
          #if (outliers.include?(-1))
          #  raise "Query is an outlier."
          #end
          temp_dm = []; temp_acts = []
          data_matrix.to_a.each_with_index { |elem, idx| temp_dm << elem unless outliers.include? idx }
          nr_cases, nr_features = temp_dm.size, temp_dm[0].size
          data_matrix = GSL::Matrix.alloc(temp_dm.flatten, nr_cases, nr_features)
          acts.each_with_index { |elem, idx| temp_acts << elem unless outliers.include? idx }
          acts = temp_acts # same nr_features


          # optimize selection of training instances -- changes cases; adjust acts accordingly
          @r.eval 'best <- vector(mode="list", length=5)'
          @r.x = data_matrix.to_a.flatten
          @r.y = acts.to_a.flatten
          @r.eval "x <- matrix(x, #{nr_cases}, #{nr_features}, byrow=T)"
          @r.eval 'df <- data.frame(y,x)'
          @r.eval 'fstr <- "y ~ ."'
          @r.eval "fit <- mvr( formula = as.formula(fstr), data=df, method = \"kernelpls\", validation = \"LOO\", ncomp=#{[ (nr_cases / 3).floor, nr_features ].min})" # was using: ncomp=#{maxcols}
          @r.eval 'rmseLoo <- matrix( RMSEP( fit, "CV" )$val )'
          @r.eval 'r2Loo <- matrix( R2( fit, "CV" )$val )'
          LOGGER.debug "RMSE (internal LOO): #{@r.rmseLoo.to_a.flatten.collect { |v| sprintf("%.2f", v) }.join(", ") }"
          @r.eval 'ncompLoo <- which( rmseLoo<=quantile(rmseLoo,.1) )[1]' # get min RMSE (10% quantile)
          @r.eval " 
            best[[1]] = #{nr_cases}
            best[[2]] = ncompLoo
            best[[3]] = rmseLoo[ncompLoo]
            best[[4]] = fit
            best[[5]] = r2Loo[ncompLoo]
          "

          # build model on best selection
          @r.eval 'best_values = c(best[[1]], best[[2]], best[[3]], best[[5]])' # Ruby-index: 0-nr neighbors, 1-nr components, 2-RMSE, 3-R2
                                                                                # Must use plain value ruby array, otherwise rinruby fails
          LOGGER.debug "Model based on #{@r.best_values[0].to_i} neighbors and #{@r.best_values[1].to_i} components, RMSE #{sprintf("%.2f", @r.best_values[2])} R2 #{sprintf("%.2f", @r.best_values[3])}."
          @r.q = query_matrix.to_a.flatten
          @r.eval "q <- data.frame( matrix( q, 1 ,#{nr_features} ) )"
          @r.eval 'names(q) = names(df)[2:length(names(df))]'
          @r.eval 'pred <- drop( predict( best[[4]], newdata = q, ncomp=best[[2]] ) )'
          @r.pred.to_a.flatten[0] # [1] lwr, [2] upr confidence limit NOT IMPLEMENTED.

        rescue Exception => e
          LOGGER.debug "#{e.class}: #{e.message}"
          #LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
        end

      end



      # Classification with majority vote from neighbors weighted by similarity
      # @param [Hash] params Keys `:acts, :sims, :value_map` are required
      # @return [Numeric] A prediction value.
      def self.weighted_majority_vote(params)

        neighbor_contribution = 0.0
        confidence_sum = 0.0
        confidence = 0.0
        prediction = nil

        LOGGER.debug "Weighted Majority Vote Classification."
        params[:acts].each_index do |idx|
          neighbor_weight = params[:sims][1][idx]
          neighbor_contribution += params[:acts][idx] * neighbor_weight

          if params[:value_map].size == 2 # AM: provide compat to binary classification: 1=>false 2=>true
            case params[:acts][idx]
            when 1
              confidence_sum -= neighbor_weight
            when 2
              confidence_sum += neighbor_weight
            end
          else
            confidence_sum += neighbor_weight
          end
        end

        if params[:value_map].size == 2 
          if confidence_sum >= 0.0
            prediction = 2 unless params[:acts].size==0
          elsif confidence_sum < 0.0
            prediction = 1 unless params[:acts].size==0
          end
        else 
          prediction = (neighbor_contribution/confidence_sum).round  unless params[:acts].size==0  # AM: new multinomial prediction
        end 
        LOGGER.debug "Prediction is: '" + prediction.to_s + "'." unless prediction.nil?
        confidence = confidence_sum/params[:acts].size if params[:acts].size > 0
        LOGGER.debug "Confidence is: '" + confidence.to_s + "'." unless prediction.nil?
        return {:prediction => prediction, :confidence => confidence.abs}
      end



      # Local support vector regression from neighbors 
      # @param [Hash] params Keys `:props, :acts, :sims` are required
      # @return [Numeric] A prediction value.
      def self.local_svm_regression(params)

        begin
          confidence = 0.0
          prediction = nil

          LOGGER.debug "Local SVM Regression."
          if params[:acts].size>0

            ## Transform data (discussion: http://goo.gl/U8Klu)
            if params[:props]
              n_prop = params[:props][0].collect
              q_prop = params[:props][1].collect
              props = [ n_prop, q_prop ]
            end

            # Transform y
            acts_autoscaler = OpenTox::Transform::LogAutoScale.new(params[:acts].to_gv)
            acts = acts_autoscaler.vs.to_a

            # Predict
            prediction = params[:props] ? local_svm_prop( props, acts, "nu-svr") : local_svm( params[:sims], acts, "nu-svr")

            # Restore
            prediction = acts_autoscaler.restore( [ prediction ].to_gv )[0]
            prediction = nil if prediction.infinite? 
            LOGGER.debug "Prediction is: '" + prediction.to_s + "'."
            params[:conf_stdev] = false if params[:conf_stdev].nil?
            confidence = get_confidence({:sims => params[:sims][1], :acts => params[:acts], :conf_stdev => params[:conf_stdev]})
            confidence = 0.0 if prediction.nil?

          end

          {:prediction => prediction, :confidence => confidence}
        rescue Exception => e
          LOGGER.debug "#{e.class}: #{e.message}"
          LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
        end
        
      end

      # Local support vector classification from neighbors 
      # @param [Hash] params Keys `:props, :acts, :sims` are required
      # @return [Numeric] A prediction value.
      def self.local_svm_classification(params)

        confidence = 0.0
        prediction = nil

        LOGGER.debug "Local SVM Classification."
        if params[:acts].size>0
          prediction = params[:props] ? local_svm_prop( params[:props], params[:acts], "C-bsvc") : local_svm( params[:sims], params[:acts], "C-bsvc") 
          LOGGER.debug "Prediction is: '" + prediction.to_s + "'."
          params[:conf_stdev] = false if params[:conf_stdev].nil?
          confidence = get_confidence({:sims => params[:sims][1], :acts => params[:acts], :conf_stdev => params[:conf_stdev]})
          confidence = 0.0 if prediction.nil?
        end
        {:prediction => prediction, :confidence => confidence}
        
      end


      # Local support vector prediction from neighbors. 
      # Uses pre-defined Kernel Matrix.
      # Not to be called directly (use local_svm_regression or local_svm_classification).
      # @param [Array] sims, similarities for neighbors.
      # @param [Array] acts, activities for neighbors.
      # @param [String] type, one of "nu-svr" (regression) or "C-bsvc" (classification).
      # @return [Numeric] A prediction value.
      def self.local_svm(sims, acts, type)
        LOGGER.debug "Local SVM (Weighted Tanimoto Kernel)."

        gram_matrix = [] # square matrix of similarities between neighbors; implements weighted tanimoto kernel
        
        prediction = nil
        if Algorithm::zero_variance? acts
          prediction = acts[0]
        else
          #LOGGER.debug gram_matrix.to_yaml
          @r = RinRuby.new(false,false) # global R instance leads to Socket errors after a large number of requests
          @r.eval "library('kernlab')" # this requires R package "kernlab" to be installed
          LOGGER.debug "Setting R data ..."
          # set data
          @r.gram_matrix = sims[0].flatten
          @r.sims = sims[1]
          @r.n = acts.size
          @r.y = acts

          begin
            LOGGER.debug "Preparing R data ..."
            # prepare data
            @r.eval "y<-as.vector(y)"
            @r.eval "gram_matrix<-as.kernelMatrix(matrix(gram_matrix,n,n))"
            @r.eval "sims<-as.vector(sims)"
            
            # model + support vectors
            LOGGER.debug "Creating SVM model ..."
            @r.eval "model<-ksvm(gram_matrix, y, kernel=matrix, type=\"#{type}\", nu=0.5)"
            @r.eval "sv<-as.vector(SVindex(model))"
            @r.eval "sims<-sims[sv]"
            @r.eval "sims<-as.kernelMatrix(matrix(sims,1))"
            LOGGER.debug "Predicting ..."
            if type == "nu-svr" 
              @r.eval "p<-predict(model,sims)[1,1]"
            elsif type == "C-bsvc"
              @r.eval "p<-predict(model,sims)"
            end
            if type == "nu-svr"
              prediction = @r.p
            elsif type == "C-bsvc"
              #prediction = (@r.p.to_f == 1.0 ? true : false)
              prediction = @r.p
            end
            @r.quit # free R
          rescue Exception => e
            LOGGER.debug "#{e.class}: #{e.message}"
            LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
          end

        end
        prediction
      end

      # Local support vector prediction from neighbors. 
      # Uses propositionalized setting.
      # Not to be called directly (use local_svm_regression or local_svm_classification).
      # @param [Array] props, propositionalization of neighbors and query structure e.g. [ Array_for_q, two-nested-Arrays_for_n ]
      # @param [Array] acts, activities for neighbors.
      # @param [String] type, one of "nu-svr" (regression) or "C-bsvc" (classification).
      # @return [Numeric] A prediction value.
      def self.local_svm_prop(props, acts, type)

          LOGGER.debug "Local SVM (Propositionalization / Kernlab Kernel)."
          n_prop = props[0] # is a matrix, i.e. two nested Arrays.
          q_prop = props[1] # is an Array.

          prediction = nil
          if Algorithm::zero_variance? acts
            prediction = acts[0]
          else
            #LOGGER.debug gram_matrix.to_yaml
            @r = RinRuby.new(false,false) # global R instance leads to Socket errors after a large number of requests
            @r.eval "library('kernlab')" # this requires R package "kernlab" to be installed
            LOGGER.debug "Setting R data ..."
            # set data
            @r.n_prop = n_prop.flatten
            @r.n_prop_x_size = n_prop.size
            @r.n_prop_y_size = n_prop[0].size
            @r.y = acts
            @r.q_prop = q_prop

            begin
              LOGGER.debug "Preparing R data ..."
              # prepare data
              @r.eval "y<-matrix(y)"
              @r.eval "prop_matrix<-matrix(n_prop, n_prop_x_size, n_prop_y_size, byrow=TRUE)"
              @r.eval "q_prop<-matrix(q_prop, 1, n_prop_y_size, byrow=TRUE)"
              
              # model + support vectors
              LOGGER.debug "Creating SVM model ..."
              @r.eval "model<-ksvm(prop_matrix, y, type=\"#{type}\", nu=0.5)"
              LOGGER.debug "Predicting ..."
              if type == "nu-svr" 
                @r.eval "p<-predict(model,q_prop)[1,1]"
              elsif type == "C-bsvc"
                @r.eval "p<-predict(model,q_prop)"
              end
              prediction = @r.p
              @r.quit # free R
            rescue Exception => e
              LOGGER.debug "#{e.class}: #{e.message}"
              LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
            end
          end
          prediction
      end

    end
    
    module Substructure
      include Algorithm
      # Substructure matching
      # @param [Hash] required keys: compound, features
      # @return [Array] Array with matching Smarts
      def self.match(params)
        params[:compound].match(params[:features])
      end
      
      # Substructure matching with number of non-unique hits
      # @param [Hash] required keys: compound, features
      # @return [Hash] Hash with matching Smarts and number of hits 
      def self.match_hits(params)
        params[:compound].match_hits(params[:features])
      end
      
      # Substructure matching with number of non-unique hits
      # @param [Hash] required keys: compound, features, feature_dataset_uri, pc_type
      # @return [Hash] Hash with matching Smarts and number of hits 
      def self.lookup(params)
        params[:compound].lookup(params[:features], params[:feature_dataset_uri],params[:pc_type])
      end  
    end

    module Dataset
      include Algorithm
      # API should match Substructure.match
      def features(dataset_uri,compound_uri)
      end
    end
  end
end
