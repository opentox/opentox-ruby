# R integration
# workaround to initialize R non-interactively (former rinruby versions did this by default)
# avoids compiling R with X
R = nil
require "rinruby" 
require "statsample"
require 'uri'

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
      # @param [Hash] params Keys `:neighbors,:compound,:features,:p_values,:similarity_algorithm,:prop_kernel,:value_map,:transform` are required
      # @return [Numeric] A prediction value.
      def self.local_mlr_prop(params)

        confidence=0.0
        prediction=nil

        if params[:neighbors].size>0
          props = params[:prop_kernel] ? get_props(params) : nil
          acts = params[:neighbors].collect { |n| act = n[:activity].to_f }
          sims = params[:neighbors].collect { |n| Algorithm.gauss(n[:similarity]) }
          LOGGER.debug "Local MLR (Propositionalization / GSL)."
          prediction = mlr( {:n_prop => props[0], :q_prop => props[1], :sims => sims, :acts => acts} )
          transformer = eval("OpenTox::Algorithm::Transform::#{params[:transform]["class"]}.new ([#{prediction}], #{params[:transform]["offset"]})")
          prediction = transformer.values[0]
          prediction = nil if prediction.infinite? || params[:prediction_min_max][1] < prediction || params[:prediction_min_max][0] > prediction  
          LOGGER.debug "Prediction is: '" + prediction.to_s + "'."
          params[:conf_stdev] = false if params[:conf_stdev].nil?
          confidence = get_confidence({:sims => sims, :acts => acts, :neighbors => params[:neighbors], :conf_stdev => params[:conf_stdev]})
          confidence = nil if prediction.nil?
        end
        {:prediction => prediction, :confidence => confidence}

      end

      # Multi-linear regression weighted by similarity.
      # Objective Feature Selection, Principal Components Analysis, Scaling of Axes.
      # @param [Hash] params Keys `:n_prop, :q_prop, :sims, :acts` are required
      # @return [Numeric] A prediction value.
      def self.mlr(params)

        # GSL matrix operations: 
        # to_a : row-wise conversion to nested array
        #
        # Statsample operations (build on GSL):
        # to_scale: convert into Statsample format

        begin
          n_prop = params[:n_prop].collect { |v| v }
          q_prop = params[:q_prop].collect { |v| v }
          n_prop << q_prop # attach q_prop
          nr_cases, nr_features = get_sizes n_prop
          data_matrix = GSL::Matrix.alloc(n_prop.flatten, nr_cases, nr_features)

          # Principal Components Analysis
          LOGGER.debug "PCA..."
          pca = OpenTox::Algorithm::Transform::PCA.new(data_matrix)
          data_matrix = pca.data_transformed_matrix

          # Attach intercept column to data
          intercept = GSL::Matrix.alloc(Array.new(nr_cases,1.0),nr_cases,1)
          data_matrix = data_matrix.horzcat(intercept)
          (0..data_matrix.size2-2).each { |i|
            autoscaler = OpenTox::Algorithm::Transform::AutoScale.new(data_matrix.col(i))
            data_matrix.col(i)[0..data_matrix.size1-1] = autoscaler.scaled_values
          }

          # Detach query instance
          n_prop = data_matrix.to_a
          q_prop = n_prop.pop 
          nr_cases, nr_features = get_sizes n_prop
          data_matrix = GSL::Matrix.alloc(n_prop.flatten, nr_cases, nr_features)

          # model + support vectors
          LOGGER.debug "Creating MLR model ..."
          c, cov, chisq, status = GSL::MultiFit::wlinear(data_matrix, params[:sims].to_scale.to_gsl, params[:acts].to_scale.to_gsl)
          GSL::MultiFit::linear_est(q_prop.to_scale.to_gsl, c, cov)[0]
        rescue Exception => e
          LOGGER.debug "#{e.class}: #{e.message}"
        end

      end

      # Classification with majority vote from neighbors weighted by similarity
      # @param [Hash] params Keys `:neighbors,:compound,:features,:p_values,:similarity_algorithm,:prop_kernel,:value_map,:transform` are required
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
      # @param [Hash] params Keys `:neighbors,:compound,:features,:p_values,:similarity_algorithm,:prop_kernel,:value_map,:transform` are required
      # @return [Numeric] A prediction value.
      def self.local_svm_regression(params)

        confidence = 0.0
        prediction = nil
        if params[:neighbors].size>0
          props = params[:prop_kernel] ? get_props(params) : nil
          acts = params[:neighbors].collect{ |n| n[:activity].to_f }
          sims = params[:neighbors].collect{ |n| Algorithm.gauss(n[:similarity]) }
          prediction = props.nil? ? local_svm(acts, sims, "nu-svr", params) : local_svm_prop(props, acts, "nu-svr")
          transformer = eval("OpenTox::Algorithm::Transform::#{params[:transform]["class"]}.new ([#{prediction}], #{params[:transform]["offset"]})")
          prediction = transformer.values[0]
          prediction = nil if prediction.infinite? || params[:prediction_min_max][1] < prediction || params[:prediction_min_max][0] > prediction  
          LOGGER.debug "Prediction is: '" + prediction.to_s + "'."
          params[:conf_stdev] = false if params[:conf_stdev].nil?
          confidence = get_confidence({:sims => sims, :acts => acts, :neighbors => params[:neighbors], :conf_stdev => params[:conf_stdev]})
          confidence = nil if prediction.nil?
        end
        {:prediction => prediction, :confidence => confidence}
        
      end

      # Local support vector classification from neighbors 
      # @param [Hash] params Keys `:neighbors,:compound,:features,:p_values,:similarity_algorithm,:prop_kernel,:value_map,:transform` are required
      # @return [Numeric] A prediction value.
      def self.local_svm_classification(params)

        confidence = 0.0
        prediction = nil
        if params[:neighbors].size>0
          props = params[:prop_kernel] ? get_props(params) : nil
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
      # @param [Hash] params Keys `:neighbors,:compound,:features,:p_values,:similarity_algorithm,:prop_kernel,:value_map,:transform` are required
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

      # Get confidence for regression, with standard deviation of neighbor activity if conf_stdev is set.
      # @param[Hash] Required keys: :sims, :acts, :neighbors, :conf_stdev
      # @return[Float] Confidence
      def self.get_confidence(params)
        if params[:conf_stdev]
          sim_median = params[:sims].to_scale.median
          if sim_median.nil?
            confidence = nil
          else
            standard_deviation = params[:acts].to_scale.standard_deviation_sample
            confidence = (sim_median*Math.exp(-1*standard_deviation)).abs
            if confidence.nan?
              confidence = nil
            end
          end
        else
          conf = params[:sims].inject{|sum,x| sum + x }
          confidence = conf/params[:neighbors].size
        end
        LOGGER.debug "Confidence is: '" + confidence.to_s + "'."
        return confidence
      end

      # Get X and Y size of a nested Array (Matrix)
      def self.get_sizes(matrix)
        begin
          nr_cases = matrix.size
          nr_features = matrix[0].size
        rescue Exception => e
          LOGGER.debug "#{e.class}: #{e.message}"
          LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
        end
        #puts "NRC: #{nr_cases}, NRF: #{nr_features}"
        [ nr_cases, nr_features ]
      end

      # Calculate the propositionalization matrix aka instantiation matrix (0/1 entries for features)
      # Same for the vector describing the query compound
      # @param[Array] neighbors.
      # @param[OpenTox::Compound] query compound.
      # @param[Array] Dataset Features.
      # @param[Array] Fingerprints of neighbors.
      # @param[Float] p-values of Features.
      def self.get_props (params)
        matrix = Array.new
        begin 
          params[:neighbors].each do |n|
            n = n[:compound]
            row = []
            params[:features].each do |f|
              if ! params[:fingerprints][n].nil? 
                row << (params[:fingerprints][n].include?(f) ? (params[:p_values][f] * params[:fingerprints][n][f]) : 0.0)
              else
                row << 0.0
              end
            end
            matrix << row
          end
          row = []
          params[:features].each do |f|
            if params[:nr_hits]
              compound_feature_hits = params[:compound].match_hits([f])
              row << (compound_feature_hits.size == 0 ? 0.0 : (params[:p_values][f] * compound_feature_hits[f]))
            else
              row << (params[:compound].match([f]).size == 0 ? 0.0 : params[:p_values][f])
            end
          end
        rescue Exception => e
          LOGGER.debug "get_props failed with '" + $! + "'"
        end
        [ matrix, row ]
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

    module Transform
      include Algorithm

      # The transformer that inverts values.
      # 1/x is used, after values have been moved >= 1.
      class Inverter
        attr_accessor :offset, :values

        # @params[Array] Values to transform.
        # @params[Float] Offset for restore.
        def initialize *args
          case args.size
          when 1
            begin
              values=args[0]
              raise "Cannot transform, values empty." if @values.size==0
              @values = values.collect { |v| -1.0 * v }  
              @offset = 1.0 - @values.minmax[0] 
              @offset = -1.0 * @offset if @offset>0.0 
              @values.collect! { |v| v - @offset }   # slide >1
              @values.collect! { |v| 1 / v }         # invert to [0,1]
            rescue Exception => e
              LOGGER.debug "#{e.class}: #{e.message}"
              LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
            end
          when 2
            @offset = args[1].to_f
            @values = args[0].collect { |v| 1 / v }
            @values.collect! { |v| v + @offset }
            @values.collect! { |v| -1.0 * v }
          end
        end
      end

      # The transformer that takes logs.
      # Log10 is used, after values have been moved > 0.
      class Log10
        attr_accessor :offset, :values

        # @params[Array] Values to transform / restore.
        # @params[Float] Offset for restore.
        def initialize *args
          @distance_to_zero = 0.000000001 # 1 / 1 billion
          case args.size
          when 1
            begin
              values=args[0]
              raise "Cannot transform, values empty." if values.size==0
              @offset = values.minmax[0] 
              @offset = -1.0 * @offset if @offset>0.0 
              @values = values.collect { |v| v - @offset }   # slide > anchor
              @values.collect! { |v| v + @distance_to_zero }  #
              @values.collect! { |v| Math::log10 v } # log10 (can fail)
            rescue Exception => e
              LOGGER.debug "#{e.class}: #{e.message}"
              LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
            end
          when 2
            @offset = args[1].to_f
            @values = args[0].collect { |v| 10**v }
            @values.collect! { |v| v - @distance_to_zero }
            @values.collect! { |v| v + @offset }
          end
        end
      end

      # The transformer that does nothing (No OPeration).
      class NOP
        attr_accessor :offset, :values

        # @params[Array] Values to transform / restore.
        # @params[Float] Offset for restore.
        def initialize *args
          @offset = 0.0
          @distance_to_zero = 0.0
          case args.size
          when 1
            @values = args[0]
          when 2
            @values = args[0]
          end
        end
      end


      # Auto-Scaler for Arrays
      # Center on mean and divide by standard deviation
      class AutoScale 
        attr_accessor :scaled_values, :mean, :stdev

        # @params[Array] Values to transform.
        def initialize values
          @scaled_values = values
          @mean = @scaled_values.to_scale.mean
          @stdev = @scaled_values.to_scale.standard_deviation_sample
          @scaled_values = @scaled_values.collect {|vi| vi - @mean }
          @scaled_values.collect! {|vi| vi / @stdev } unless @stdev == 0.0
        end
      end

      # Principal Components Analysis
      # Statsample Library (http://ruby-statsample.rubyforge.org/) by C. Bustos
      class PCA
        attr_accessor :data_matrix, :data_transformed_matrix, :eigenvector_matrix, :eigenvalue_sums, :autoscaler

        # Creates a transformed dataset as GSL::Matrix.
        # @param [GSL::Matrix] Data matrix.
        # @param [Float] Compression ratio from [0,1].
        # @return [GSL::Matrix] Data transformed matrix.
        def initialize data_matrix, compression=0.05
          begin
            @data_matrix = data_matrix
            @compression = compression.to_f
            @stdev = Array.new
            @mean = Array.new

            # Objective Feature Selection
            raise "Error! PCA needs at least two dimensions." if data_matrix.size2 < 2
            @data_matrix_selected = nil
            (0..@data_matrix.size2-1).each { |i|
              if !Algorithm::zero_variance?(@data_matrix.col(i).to_a)
                if @data_matrix_selected.nil?
                  @data_matrix_selected = GSL::Matrix.alloc(@data_matrix.size1, 1) 
                  @data_matrix_selected.col(0)[0..@data_matrix.size1-1] = @data_matrix.col(i)
                else
                  @data_matrix_selected = @data_matrix_selected.horzcat(GSL::Matrix.alloc(@data_matrix.col(i).to_a,@data_matrix.size1, 1))
                end
              end             
            }
            raise "Error! PCA needs at least two dimensions." if (@data_matrix_selected.nil? || @data_matrix_selected.size2 < 2)

            # Scaling of Axes
            @data_matrix_scaled = GSL::Matrix.alloc(@data_matrix_selected.size1, @data_matrix_selected.size2)
            (0..@data_matrix_selected.size2-1).each { |i|
              @autoscaler = OpenTox::Algorithm::Transform::AutoScale.new(@data_matrix_selected.col(i))
              @data_matrix_scaled.col(i)[0..@data_matrix.size1-1] = @autoscaler.scaled_values
              @stdev << @autoscaler.stdev
              @mean << @autoscaler.mean
            }

            data_matrix_hash = Hash.new
            (0..@data_matrix_scaled.size2-1).each { |i|
              column_view = @data_matrix_scaled.col(i)
              data_matrix_hash[i] = column_view.to_scale
            }
            dataset_hash = data_matrix_hash.to_dataset # see http://goo.gl/7XcW9
            cor_matrix=Statsample::Bivariate.correlation_matrix(dataset_hash)
            pca=Statsample::Factor::PCA.new(cor_matrix)
            pca.eigenvalues.each { |ev| raise "PCA failed!" unless !ev.nan? }
            @eigenvalue_sums = Array.new
            (0..dataset_hash.fields.size-1).each { |i|
              @eigenvalue_sums << pca.eigenvalues[0..i].inject{ |sum, ev| sum + ev }
            }
            eigenvectors_selected = Array.new
            pca.eigenvectors.each_with_index { |ev, i|
              if (@eigenvalue_sums[i] <= ((1.0-@compression)*dataset_hash.fields.size)) || (eigenvectors_selected.size == 0)
                eigenvectors_selected << ev.to_a
              end
            }
            @eigenvector_matrix = GSL::Matrix.alloc(eigenvectors_selected.flatten, eigenvectors_selected.size, dataset_hash.fields.size).transpose
            dataset_matrix = dataset_hash.to_gsl.transpose
            @data_transformed_matrix = (@eigenvector_matrix.transpose * dataset_matrix).transpose
          rescue Exception => e
              LOGGER.debug "#{e.class}: #{e.message}"
              LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
          end
        end

        # Restores data in the original feature space (possibly with compression loss).
        # @return [GSL::Matrix] Data matrix.
        def restore
          begin 
            data_matrix_restored = (@eigenvector_matrix * @data_transformed_matrix.transpose).transpose # reverse pca
            # reverse scaling
            (0..data_matrix_restored.size2-1).each { |i|
              data_matrix_restored.col(i)[0..data_matrix_restored.size1-1] *= @stdev[i] unless @stdev[i] == 0.0
              data_matrix_restored.col(i)[0..data_matrix_restored.size1-1] += @mean[i]
            }
            data_matrix_restored
          rescue Exception => e
            LOGGER.debug "#{e.class}: #{e.message}"
            LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
          end
        end

      end

    end
    
    # Gauss kernel
    # @return [Float] 
    def self.gauss(x, sigma = 0.3) 
      d = 1.0 - x.to_f
      Math.exp(-(d*d)/(2*sigma*sigma))
    end

    # For symbolic features
    # @param [Array] Array to test, must indicate non-occurrence with 0.
    # @return [Boolean] Whether the feature is singular or non-occurring or present everywhere.
    def self.isnull_or_singular?(array)
      nr_zeroes = array.count(0)
      return (nr_zeroes == array.size) ||    # remove non-occurring feature
             (nr_zeroes == array.size-1) ||  # remove singular feature
             (nr_zeroes == 0)                # also remove feature present everywhere
    end

    # Numeric value test
    # @param[Object] value
    # @return [Boolean] Whether value is a number
    def self.numeric?(value)
      true if Float(value) rescue false
    end

    # For symbolic features
    # @param [Array] Array to test, must indicate non-occurrence with 0.
    # @return [Boolean] Whether the feature has variance zero.
    def self.zero_variance?(array)
      return (array.to_scale.variance_sample == 0.0)
    end
    
    # Sum of an array for Arrays.
    # @param [Array] Array with values
    # @return [Integer] Sum of size of values
    def self.sum_size(array)
      sum=0
      array.each { |e| sum += e.size }
      return sum
    end

    # Minimum Frequency
    # @param [Integer] per-mil value
    # return [Integer] min-frequency
    def self.min_frequency(training_dataset,per_mil)
      minfreq = per_mil * training_dataset.compounds.size.to_f / 1000.0 # AM sugg. 8-10 per mil for BBRC, 50 per mil for LAST
      minfreq = 2 unless minfreq > 2
      Integer (minfreq)
    end

    # Effect calculation for classification
    # @param [Array] Array of occurrences per class in the form of Enumerables.
    # @param [Array] Array of database instance counts per class.
    def self.effect(occurrences, db_instances)
      max=0
      max_value=0
      nr_o = self.sum_size(occurrences)
      nr_db = db_instances.to_scale.sum

      occurrences.each_with_index { |o,i| # fminer outputs occurrences sorted reverse by activity.
        actual = o.size.to_f/nr_o
        expected = db_instances[i].to_f/nr_db
        if actual > expected
          if ((actual - expected) / actual) > max_value
           max_value = (actual - expected) / actual # 'Schleppzeiger'
            max = i
          end
        end
      }
      max
    end
    
    # Returns Support value of an fingerprint
    # @param [Hash] params Keys: `:compound_features_hits, :weights, :training_compound_features_hits, :features, :nr_hits:, :mode` are required
    # return [Numeric] Support value 
    def self.p_sum_support(params)
      p_sum = 0.0
        params[:features].each{|f|
        compound_hits = params[:compound_features_hits][f]
        neighbor_hits = params[:training_compound_features_hits][f] 
        p_sum += eval("(Algorithm.gauss(params[:weights][f]) * ([compound_hits, neighbor_hits].compact.#{params[:mode]}))")
      }
      p_sum 
    end
                
  end
end


