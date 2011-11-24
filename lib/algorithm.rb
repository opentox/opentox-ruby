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

    # Utility methods without dedicated webservices

    # Similarity calculations
    module Similarity
      include Algorithm

      # Tanimoto similarity
      # @param [Array] features_a Features of first compound
      # @param [Array] features_b Features of second compound
      # @param [optional, Hash] weights Weights for all features
      # @param [optional, Hash] params Keys: `:training_compound, :compound, :training_compound_features_hits, :nr_hits, :compound_features_hits` are required
      # @return [Float] (Weighted) tanimoto similarity
      def self.tanimoto(features_a,features_b,weights=nil,params=nil)
        common_features = features_a & features_b
        all_features = (features_a + features_b).uniq
        #LOGGER.debug "dv --------------- common: #{common_features}, all: #{all_features}"
        if common_features.size > 0
          if weights
            #LOGGER.debug "nr_hits: #{params[:nr_hits]}"
            if !params.nil? && params[:nr_hits]
              params[:weights] = weights
              params[:mode] = "min"
              params[:features] = common_features
              common_p_sum = Algorithm.p_sum_support(params)
              params[:mode] = "max"
              params[:features] = all_features
              all_p_sum = Algorithm.p_sum_support(params)
            else
              common_p_sum = 0.0
              common_features.each{|f| common_p_sum += Algorithm.gauss(weights[f])}
              all_p_sum = 0.0
              all_features.each{|f| all_p_sum += Algorithm.gauss(weights[f])}
            end
            #LOGGER.debug "common_p_sum: #{common_p_sum}, all_p_sum: #{all_p_sum}, c/a: #{common_p_sum/all_p_sum}"
            common_p_sum/all_p_sum
          else
            #LOGGER.debug "common_features : #{common_features}, all_features: #{all_features}, c/a: #{(common_features.size/all_features.size).to_f}"
            common_features.size.to_f/all_features.size.to_f
          end
        else
          0.0
        end
      end

      # Euclidean similarity
      # @param [Hash] properties_a Properties of first compound
      # @param [Hash] properties_b Properties of second compound
      # @param [optional, Hash] weights Weights for all properties
      # @return [Float] (Weighted) euclidean similarity
      def self.euclidean(properties_a,properties_b,weights=nil)
        common_properties = properties_a.keys & properties_b.keys
        if common_properties.size > 1
          dist_sum = 0
          common_properties.each do |p|
            if weights
              dist_sum += ( (properties_a[p] - properties_b[p]) * Algorithm.gauss(weights[p]) )**2
            else
              dist_sum += (properties_a[p] - properties_b[p])**2
            end
          end
          1/(1+Math.sqrt(dist_sum))
        else
          0.0
        end
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
      # @param [Hash] params Keys `:neighbors,:compound,:features,:p_values,:similarity_algorithm,:prop_kernel,:value_map` are required.
      # @return [Numeric] A prediction value.

      def self.local_mlr_prop(params)

        confidence=0.0
        prediction=nil

        if params[:neighbors].size>0
          props, ids = params[:prop_kernel] ? get_props_pc(params) : nil
          acts = params[:neighbors].collect { |n| n[:activity].to_f }
          acts = acts.collect { |e| e if ids.include? acts.index(e) }.compact # remove acts of removed neighbors
          sims = params[:neighbors].collect { |n| Algorithm.gauss(n[:similarity]) }
          sims = sims.collect { |e| e if ids.include? sims.index(e) }.compact # remove acts of removed neighbors
          maxcols = ( params[:maxcols].nil? ? (sims.size/3.0).ceil : params[:maxcols] )
          LOGGER.debug "Local MLR (Propositionalization / GSL)."
          prediction = pcr( {:n_prop => props[0], :q_prop => props[1], :sims => sims, :acts => acts, :maxcols => maxcols} )
          prediction = nil if prediction.infinite? || params[:prediction_min_max][1] < prediction || params[:prediction_min_max][0] > prediction  
          LOGGER.debug "Prediction is: '" + prediction.to_s + "'."
          params[:conf_stdev] = false if params[:conf_stdev].nil?
          confidence = get_confidence({:sims => sims, :acts => acts, :neighbors => params[:neighbors], :conf_stdev => params[:conf_stdev]})
          confidence = nil if prediction.nil?
        end
        {:prediction => prediction, :confidence => confidence}

      end

      # Multi-linear regression (unweighted).
      # Objective Feature Selection, Scaling of Axes, Principal Components Analysis.
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
          maxcols = params[:maxcols]

          nr_cases, nr_features = get_sizes n_prop
          data_matrix = GSL::Matrix.alloc(n_prop.flatten, nr_cases, nr_features)
          query_matrix = GSL::Matrix.alloc(q_prop.flatten, 1, nr_features) # same nr_features

          ### Transform data (discussion: http://goo.gl/U8Klu)
          # Standardize data (scale and center), adjust query accordingly
          LOGGER.debug "Standardize..."
          (0..nr_features-1).each { |i|
            autoscaler = OpenTox::Transform::LogAutoScale.new(data_matrix.col(i))
            data_matrix.col(i)[0..nr_cases-1] = autoscaler.vs
            query_matrix.col(i)[0] = autoscaler.transform(query_matrix.col(i))[0]
          }

          # Rotate data (pca), adjust query accordingly
          LOGGER.debug "PCA..."
          pca = OpenTox::Transform::PCA.new(data_matrix,0.05,maxcols)
          data_matrix = pca.data_transformed_matrix
          nr_cases, nr_features = get_sizes data_matrix.to_a
          query_matrix = pca.transform(query_matrix)

          # Transform y
          acts_autoscaler = OpenTox::Transform::LogAutoScale.new(acts.to_gv)
          acts = acts_autoscaler.vs.to_a
          ### End of transform
          
          ### Model
          LOGGER.debug "MLR..."
          @r = RinRuby.new(false,false)   # global R instance leads to Socket errors after a large number of requests
          @r.x = data_matrix.to_a.flatten
          @r.q = query_matrix.to_a.flatten
          @r.y = acts.to_a.flatten

          @r.eval "x <- matrix(x, #{nr_cases}, #{nr_features}, byrow=T)"
          @r.eval "df <- data.frame(y,x)"

          @r.eval "fstr <- \"y ~ .\""
          @r.eval "fit <- lm( as.formula(fstr), data=df)"

          @r.eval "q <- data.frame( matrix( q, 1 ,#{nr_features} ) )"
          @r.eval "names(q) = names(df)[2:length(names(df))]"
          @r.eval "pred <- predict(fit, q, interval=\"confidence\")"
          ### End of Model
          
          point_prediction = [ @r.pred ].to_a.flatten[0] # [1] is lwr, [2] upr confidence limit.
          acts_autoscaler.restore( [ point_prediction ].to_gv )[0] # return restored value of type numeric

        rescue Exception => e
          LOGGER.debug "#{e.class}: #{e.message}"
          LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
        end

      end


      # Principal Components Regression (unweighted).
      # Objective Feature Selection, Scaling of Axes.
      # @param [Hash] params Keys `:n_prop, :q_prop, :sims, :acts, :maxcols` are required.
      # @return [Numeric] A prediction value.
      def self.pcr(params)

        # Uses Statsample Library (http://ruby-statsample.rubyforge.org/) by C. Bustos
        # Statsample operations build on GSL and offer an R-like access to data
        LOGGER.debug "MLR..."
        begin
          n_prop = params[:n_prop].collect
          q_prop = params[:q_prop].collect
          acts = params[:acts].collect
          maxcols = params[:maxcols]

          nr_cases, nr_features = get_sizes n_prop
          maxcols = nr_features if maxcols > nr_features

          data_matrix = GSL::Matrix.alloc(n_prop.flatten, nr_cases, nr_features)
          query_matrix = GSL::Matrix.alloc(q_prop.flatten, 1, nr_features) # same nr_features

          ### Transform data (discussion: http://goo.gl/U8Klu)
          # Standardize data (scale and center), adjust query accordingly
          LOGGER.debug "Standardize..."
          (0..nr_features-1).each { |i|
            autoscaler = OpenTox::Transform::LogAutoScale.new(data_matrix.col(i))
            data_matrix.col(i)[0..nr_cases-1] = autoscaler.vs
            query_matrix.col(i)[0] = autoscaler.transform(query_matrix.col(i))[0]
          }
          # Transform y
          acts_autoscaler = OpenTox::Transform::LogAutoScale.new(acts.to_gv)
          acts = acts_autoscaler.vs.to_a
          ### End of transform
          
          ### Model
          LOGGER.debug "PCR / PLSR..."
          @r = RinRuby.new(false,false)   # global R instance leads to Socket errors after a large number of requests
          @r.x = data_matrix.to_a.flatten
          @r.q = query_matrix.to_a.flatten
          @r.y = acts.to_a.flatten

          @r.eval "library(\"pls\")"
          @r.eval "x <- matrix(x, #{nr_cases}, #{nr_features}, byrow=T)"
          @r.eval "df <- data.frame(y,x)"

          @r.eval "fstr <- \"y ~ .\""
          @r.eval "fit <- mvr( formula = as.formula(fstr), data=df, ncomp=#{maxcols}, method = \"kernelpls\", validation = \"LOO\" )"

          @r.eval "r2Loo <- matrix( R2( fit, \"CV\" )$val )"
          LOGGER.debug "R-Squared (internal LOO CV), using 0 to #{maxcols} PCs: #{@r.r2Loo.to_a.flatten.collect { |v| sprintf("%.2f", v) }.join(", ") }"
          @r.eval "ncompLoo <- which.max(r2Loo) - 1"
          @r.eval "ncompLoo <- ncompLoo + 1" if (@r.ncompLoo.to_i == 0)

          LOGGER.debug "Best position: #{@r.ncompLoo.to_i}"

          @r.eval "q <- data.frame( matrix( q, 1 ,#{nr_features} ) )"
          @r.eval "names(q) = names(df)[2:length(names(df))]"
          @r.eval "pred <- drop( predict( fit, newdata = q, ncomp=ncompLoo ) )"
          ### End of Model
          
          point_prediction = [ @r.pred ].to_a.flatten[0] # [1] lwr, [2] upr confidence limit NOT IMPLEMENTED.
          acts_autoscaler.restore( [ point_prediction ].to_gv )[0] # return restored value of type numeric

        rescue Exception => e
          LOGGER.debug "#{e.class}: #{e.message}"
          LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
        end

      end





      # Classification with majority vote from neighbors weighted by similarity
      # @param [Hash] params Keys `:neighbors,:compound,:features,:p_values,:similarity_algorithm,:prop_kernel,:value_map` are required
      # @return [Numeric] A prediction value.
      def self.weighted_majority_vote(params)

        neighbor_contribution = 0.0
        confidence_sum = 0.0
        confidence = 0.0
        prediction = nil

        params[:neighbors].each do |neighbor|
          neighbor_weight = Algorithm.gauss(neighbor[:similarity]).to_f
          neighbor_contribution += neighbor[:activity].to_f * neighbor_weight

          if params[:value_map].size == 2 # AM: provide compat to binary classification: 1=>false 2=>true
            case neighbor[:activity]
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
            prediction = 2 unless params[:neighbors].size==0
          elsif confidence_sum < 0.0
            prediction = 1 unless params[:neighbors].size==0
          end
        else 
          prediction = (neighbor_contribution/confidence_sum).round  unless params[:neighbors].size==0  # AM: new multinomial prediction
        end 
        LOGGER.debug "Prediction is: '" + prediction.to_s + "'." unless prediction.nil?
        confidence = confidence_sum/params[:neighbors].size if params[:neighbors].size > 0
        LOGGER.debug "Confidence is: '" + confidence.to_s + "'." unless prediction.nil?
        return {:prediction => prediction, :confidence => confidence.abs}
      end

      # Local support vector regression from neighbors 
      # @param [Hash] params Keys `:neighbors,:compound,:features,:p_values,:similarity_algorithm,:prop_kernel,:value_map` are required
      # @return [Numeric] A prediction value.
      def self.local_svm_regression(params)

        confidence = 0.0
        prediction = nil
        if params[:neighbors].size>0
          props = params[:prop_kernel] ? get_props_fingerprints(params) : nil
          acts = params[:neighbors].collect{ |n| n[:activity].to_f }

          # Transform y
          acts_autoscaler = OpenTox::Transform::LogAutoScale.new(acts.to_gv)
          acts = acts_autoscaler.vs.to_a
 
          sims = params[:neighbors].collect{ |n| Algorithm.gauss(n[:similarity]) }
          prediction = props.nil? ? local_svm(acts, sims, "nu-svr", params) : local_svm_prop(props, acts, "nu-svr")
          prediction = nil if prediction.infinite? || params[:prediction_min_max][1] < prediction || params[:prediction_min_max][0] > prediction  

          acts_autoscaler.restore( [ prediction ].to_gv )[0] # return restored value of type numeric

          LOGGER.debug "Prediction is: '" + prediction.to_s + "'."
          params[:conf_stdev] = false if params[:conf_stdev].nil?
          confidence = get_confidence({:sims => sims, :acts => acts, :neighbors => params[:neighbors], :conf_stdev => params[:conf_stdev]})
          confidence = nil if prediction.nil?
        end
        {:prediction => prediction, :confidence => confidence}
        
      end

      # Local support vector classification from neighbors 
      # @param [Hash] params Keys `:neighbors,:compound,:features,:p_values,:similarity_algorithm,:prop_kernel,:value_map` are required
      # @return [Numeric] A prediction value.
      def self.local_svm_classification(params)

        confidence = 0.0
        prediction = nil
        if params[:neighbors].size>0
          props = params[:prop_kernel] ? get_props_fingerprints(params) : nil
          acts = params[:neighbors].collect { |n| act = n[:activity] }
          sims = params[:neighbors].collect{ |n| Algorithm.gauss(n[:similarity]) } # similarity values btwn q and nbors
          prediction = props.nil? ? local_svm(acts, sims, "C-bsvc", params) : local_svm_prop(props, acts, "C-bsvc")
          LOGGER.debug "Prediction is: '" + prediction.to_s + "'."
          params[:conf_stdev] = false if params[:conf_stdev].nil?
          confidence = get_confidence({:sims => sims, :acts => acts, :neighbors => params[:neighbors], :conf_stdev => params[:conf_stdev]})
        end
        {:prediction => prediction, :confidence => confidence}
        
      end


      # Local support vector prediction from neighbors. 
      # Uses pre-defined Kernel Matrix.
      # Not to be called directly (use local_svm_regression or local_svm_classification).
      # @param [Array] acts, activities for neighbors.
      # @param [Array] sims, similarities for neighbors.
      # @param [String] type, one of "nu-svr" (regression) or "C-bsvc" (classification).
      # @param [Hash] params Keys `:neighbors,:compound,:features,:p_values,:similarity_algorithm,:prop_kernel,:value_map` are required
      # @return [Numeric] A prediction value.
      def self.local_svm(acts, sims, type, params)
        LOGGER.debug "Local SVM (Weighted Tanimoto Kernel)."
        neighbor_matches = params[:neighbors].collect{ |n| n[:features] } # URIs of matches
        gram_matrix = [] # square matrix of similarities between neighbors; implements weighted tanimoto kernel

        prediction = nil
        if Algorithm::zero_variance? acts
          prediction = acts[0]
        else
          # gram matrix
          (0..(neighbor_matches.length-1)).each do |i|
            neighbor_i_hits = params[:fingerprints][params[:neighbors][i][:compound]]
            gram_matrix[i] = [] unless gram_matrix[i]
            # upper triangle
            ((i+1)..(neighbor_matches.length-1)).each do |j|
              neighbor_j_hits= params[:fingerprints][params[:neighbors][j][:compound]]
              sim_params = {}
              if params[:nr_hits]
                sim_params[:nr_hits] = true
                sim_params[:compound_features_hits] = neighbor_i_hits
                sim_params[:training_compound_features_hits] = neighbor_j_hits
              end
              sim = eval("#{params[:similarity_algorithm]}(neighbor_matches[i], neighbor_matches[j], params[:p_values], sim_params)")
              gram_matrix[i][j] = Algorithm.gauss(sim)
              gram_matrix[j] = [] unless gram_matrix[j] 
              gram_matrix[j][i] = gram_matrix[i][j] # lower triangle
            end
            gram_matrix[i][i] = 1.0
          end


          #LOGGER.debug gram_matrix.to_yaml
          @r = RinRuby.new(false,false) # global R instance leads to Socket errors after a large number of requests
          @r.eval "library('kernlab')" # this requires R package "kernlab" to be installed
          LOGGER.debug "Setting R data ..."
          # set data
          @r.gram_matrix = gram_matrix.flatten
          @r.n = neighbor_matches.size
          @r.y = acts
          @r.sims = sims

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

    end

    module Substructure
      include Algorithm
      # Substructure matching
      # @param [OpenTox::Compound] compound Compound
      # @param [Array] features Array with Smarts strings
      # @return [Array] Array with matching Smarts
      def self.match(compound,features)
        compound.match(features)
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
