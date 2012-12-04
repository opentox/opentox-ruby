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

    # Execute algorithm with parameters, consult OpenTox API and webservice documentation for acceptable parameters
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
      attr_accessor :prediction_feature, :training_dataset, :minfreq, :compounds, :db_class_sizes, :all_activities, :smi, :weight_feature

      def check_params(params,per_mil,subjectid=nil)
        raise OpenTox::NotFoundError.new "Please submit a dataset_uri." unless params[:dataset_uri] and  !params[:dataset_uri].nil?
        @training_dataset = OpenTox::Dataset.find "#{params[:dataset_uri]}", subjectid

        unless params[:prediction_feature] # try to read prediction_feature from dataset
          raise OpenTox::NotFoundError.new "Please provide a prediction_feature parameter" unless @training_dataset.features.size == 1
          prediction_feature = OpenTox::Feature.find(@training_dataset.features.keys.first,@subjectid)
          params[:prediction_feature] = prediction_feature.uri
        end
        @prediction_feature = OpenTox::Feature.find params[:prediction_feature], subjectid

        raise OpenTox::NotFoundError.new "No feature #{params[:prediction_feature]} in dataset #{params[:dataset_uri]}" unless @training_dataset.features and @training_dataset.features.include?(params[:prediction_feature])

        unless params[:weight_feature].nil?
          @weight_feature = OpenTox::Feature.find params[:weight_feature], subjectid
          raise OpenTox::NotFoundError.new "No feature #{params[:weight_feature]} in dataset #{params[:dataset_uri]}" unless @training_dataset.features and @training_dataset.features.include?(params[:weight_feature])
        end

        unless params[:min_frequency].nil? 
          # check for percentage
          if params[:min_frequency].include? "pc"
            per_mil=params[:min_frequency].gsub(/pc/,"")
            if OpenTox::Algorithm.numeric? per_mil
              per_mil = per_mil.to_i * 10
            else
              bad_request=true
            end
          # check for per-mil
          elsif params[:min_frequency].include? "pm"
            per_mil=params[:min_frequency].gsub(/pm/,"")
            if OpenTox::Algorithm.numeric? per_mil
              per_mil = per_mil.to_i
            else
              bad_request=true
            end
          # set minfreq directly
          else
            if OpenTox::Algorithm.numeric? params[:min_frequency]
              @minfreq=params[:min_frequency].to_i
              LOGGER.debug "min_frequency #{@minfreq}"
            else
              bad_request=true
            end
          end
          raise OpenTox::BadRequestError.new "Minimum frequency must be integer [n], or a percentage [n]pc, or a per-mil [n]pm , with n greater 0" if bad_request
        end
        if @minfreq.nil?
          @minfreq=OpenTox::Algorithm.min_frequency(@training_dataset,@prediction_feature,per_mil)
          LOGGER.debug "min_frequency #{@minfreq} (input was #{per_mil} per-mil)"
        end
      end

      def add_fminer_data(fminer_instance, value_map)

        id = 1 # fminer start id is not 0
        which_row=@training_dataset.compounds.inject({}) {|h,c| h[c]=0; h}

        @training_dataset.compounds.each do |compound|
          entry=@training_dataset.data_entries[compound]
          begin
            smiles = OpenTox::Compound.new(compound).to_smiles
          rescue
            LOGGER.warn "No resource for #{compound.to_s}"
            next
          end
          if smiles == '' or smiles.nil?
            LOGGER.warn "Cannot find smiles for #{compound.to_s}."
            next
          end

          entry && entry.each do |feature,values|
            if feature == @prediction_feature.uri
              value=values[which_row[compound]]
              if value.nil? 
                LOGGER.warn "No #{feature} activity for #{compound.to_s}."
              else
                if @prediction_feature.feature_type == "classification"
                  activity= value_map.invert[value].to_i # activities are mapped to 1..n
                  raise "activity should be mapped to 1..n for id '#{id}' with value '#{value}', value_map: #{value_map.inspect}" if activity==0
                  @db_class_sizes[activity-1].nil? ? @db_class_sizes[activity-1]=1 : @db_class_sizes[activity-1]+=1 # AM effect
                elsif @prediction_feature.feature_type == "regression"
                  activity= value.to_f 
                end
                begin
                  fminer_instance.AddCompound(smiles,id) if fminer_instance
                  fminer_instance.AddActivity(activity, id) if fminer_instance 
                  @all_activities[id]=activity # DV: insert global information
                  @compounds[id] = compound
                  @smi[id] = smiles
                  if ((not fminer_instance.nil?) and (not @weight_feature.nil?) and (@prediction_feature.feature_type == "classification"))
                    weight=entry[@weight_feature.uri][which_row[compound]].to_f # nil.to_f = 0
                    raise "weights should be positive for id '#{id}' with weight '#{weight}'" unless weight>0.0
                    fminer_instance.AddWeight(weight, id)
                  end
                  id += 1
                rescue Exception => e
                  LOGGER.warn "Could not add " + smiles + "\t" + value.to_s + " to fminer"
                  LOGGER.warn e.backtrace
                end
              end
              which_row[compound] += 1
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
        if !self.numeric? training_threshold || training_threshold <0 || training_threshold >1
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

        if !self.numeric? query_threshold || query_threshold <0 || query_threshold >1
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
        confidence = (confidence_sum/params[:acts].size).abs if params[:acts].size > 0
        LOGGER.debug "Confidence is: '" + confidence.to_s + "'." unless prediction.nil?
        return {:prediction => prediction, :confidence => confidence.abs}
      end



      # Local support vector regression from neighbors 
      # @param [Hash] params Keys `:props, :acts, :sims, :min_train_performance` are required
      # @return [Numeric] A prediction value.
      def self.local_svm_regression(params)

        begin
          confidence = 0.0
          prediction = nil

          LOGGER.debug "Local SVM."
          if params[:acts].size>0
            if params[:props]
              n_prop = params[:props][0].collect
              q_prop = params[:props][1].collect
              props = [ n_prop, q_prop ]
            end
            acts = params[:acts].collect
            prediction = local_svm_prop( props, acts, params[:min_train_performance]) # params[:props].nil? signals non-prop setting
            prediction = nil if (!prediction.nil? && prediction.infinite?)
            LOGGER.debug "Prediction is: '" + prediction.to_s + "'."
            confidence = get_confidence({:sims => params[:sims][1], :acts => params[:acts]})
            confidence = 0.0 if prediction.nil?
          end
          {:prediction => prediction, :confidence => confidence}
        rescue Exception => e
          LOGGER.debug "#{e.class}: #{e.message}"
          LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
        end

      end


      # Local support vector regression from neighbors 
      # @param [Hash] params Keys `:props, :acts, :sims, :min_train_performance` are required
      # @return [Numeric] A prediction value.
      def self.local_svm_classification(params)

        begin
          confidence = 0.0
          prediction = nil

          LOGGER.debug "Local SVM."
          if params[:acts].size>0
            if params[:props]
              n_prop = params[:props][0].collect
              q_prop = params[:props][1].collect
              props = [ n_prop, q_prop ]
            end
            acts = params[:acts].collect
            acts = acts.collect{|v| "Val" + v.to_s} # Convert to string for R to recognize classification
            prediction = local_svm_prop( props, acts, params[:min_train_performance]) # params[:props].nil? signals non-prop setting
            prediction = prediction.sub(/Val/,"") if prediction # Convert back to Float
            confidence = 0.0 if prediction.nil?
            LOGGER.debug "Prediction is: '" + prediction.to_s + "'."
            confidence = get_confidence({:sims => params[:sims][1], :acts => params[:acts]})
          end
          {:prediction => prediction, :confidence => confidence}
        rescue Exception => e
          LOGGER.debug "#{e.class}: #{e.message}"
          LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
        end

      end



      # Local support vector prediction from neighbors. 
      # Uses propositionalized setting.
      # Not to be called directly (use local_svm_regression or local_svm_classification).
      # @param [Array] props, propositionalization of neighbors and query structure e.g. [ Array_for_q, two-nested-Arrays_for_n ]
      # @param [Array] acts, activities for neighbors.
      # @param [Float] min_train_performance, parameter to control censoring
      # @return [Numeric] A prediction value.
      def self.local_svm_prop(props, acts, min_train_performance)

        LOGGER.debug "Local SVM (Propositionalization / Kernlab Kernel)."
        n_prop = props[0] # is a matrix, i.e. two nested Arrays.
        q_prop = props[1] # is an Array.

        prediction = nil
        if Algorithm::zero_variance? acts
          prediction = acts[0]
        else
          #LOGGER.debug gram_matrix.to_yaml
          @r = RinRuby.new(true,false) # global R instance leads to Socket errors after a large number of requests
          @r.eval "suppressPackageStartupMessages(library('caret'))" # requires R packages "caret" and "kernlab"
          @r.eval "suppressPackageStartupMessages(library('doMC'))" # requires R packages "multicore"
          @r.eval "registerDoMC()" # switch on parallel processing
          @r.eval "set.seed(1)"
          begin

            # set data
            LOGGER.debug "Setting R data ..."
            @r.n_prop = n_prop.flatten
            @r.n_prop_x_size = n_prop.size
            @r.n_prop_y_size = n_prop[0].size
            @r.y = acts
            @r.q_prop = q_prop
            #@r.eval "y = matrix(y)"
            @r.eval "prop_matrix = matrix(n_prop, n_prop_x_size, n_prop_y_size, byrow=T)"
            @r.eval "q_prop = matrix(q_prop, 1, n_prop_y_size, byrow=T)"

            # prepare data
            LOGGER.debug "Preparing R data ..."
            @r.eval <<-EOR
              weights=NULL
              if (!(class(y) == 'numeric')) { 
                y = factor(y)
                suppressPackageStartupMessages(library('class')) 
                weights=unlist(as.list(prop.table(table(y))))
                weights=(weights-1)^2
              }
            EOR

            @r.eval <<-EOR
              rem = nearZeroVar(prop_matrix)
              if (length(rem) > 0) {
                prop_matrix = prop_matrix[,-rem,drop=F]
                q_prop = q_prop[,-rem,drop=F]
              }
              rem = findCorrelation(cor(prop_matrix))
              if (length(rem) > 0) {
                prop_matrix = prop_matrix[,-rem,drop=F]
                q_prop = q_prop[,-rem,drop=F]
              }
            EOR

            # model + support vectors
            LOGGER.debug "Creating R SVM model ..."
            train_success = @r.eval <<-EOR
              model = train(prop_matrix,y,
                             method="svmradial",
                             preProcess=c("center", "scale"),
                             class.weights=weights,
                             trControl=trainControl(method="LGOCV",number=10),
                             tuneLength=8
                           )
              perf = ifelse ( class(y)!='numeric', max(model$results$Accuracy), model$results[which.min(model$results$RMSE),]$Rsquared )
            EOR


            # prediction
            LOGGER.debug "Predicting ..."
            @r.eval "p = predict(model,q_prop)"
            @r.eval "if (class(y)!='numeric') p = as.character(p)"
            prediction = @r.p

            # censoring
            prediction = nil if ( @r.perf.nan? || @r.perf < min_train_performance )
            prediction = nil unless train_success
            LOGGER.debug "Performance: #{sprintf("%.2f", @r.perf)}"
          rescue Exception => e
            LOGGER.debug "#{e.class}: #{e.message}"
            LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
          end
          @r.quit # free R
        end
        prediction
      end

    end

    module FeatureSelection
      include Algorithm
      # Recursive Feature Elimination using caret
      # @param [Hash] required keys: ds_csv_file, prediction_feature, fds_csv_file (dataset CSV file, prediction feature column name, and feature dataset CSV file), optional: del_missing (delete rows with missing values).
      # @return [String] feature dataset CSV file composed of selected features.
      def self.rfe(params)
        @r=RinRuby.new(false,false)
        @r.ds_csv_file = params[:ds_csv_file].to_s
        @r.prediction_feature = params[:prediction_feature].to_s
        @r.fds_csv_file = params[:fds_csv_file].to_s
        @r.del_missing = params[:del_missing] == true ? 1 : 0
        r_result_file = params[:fds_csv_file].sub("rfe_", "rfe_R_")
        @r.f_fds_r = r_result_file.to_s

        # need packs 'randomForest', 'RANN'
        @r.eval <<-EOR
          suppressPackageStartupMessages(library('caret'))
          suppressPackageStartupMessages(library('randomForest'))
          suppressPackageStartupMessages(library('RANN'))
          suppressPackageStartupMessages(library('doMC'))
          registerDoMC()
          set.seed(1)

          acts = read.csv(ds_csv_file, check.names=F)
          feats = read.csv(fds_csv_file, check.names=F)
          ds = merge(acts, feats, by="SMILES") # duplicates features for duplicate SMILES :-)

          features = ds[,(dim(acts)[2]+1):(dim(ds)[2])]
          y = ds[,which(names(ds) == prediction_feature)] 

          # assumes a data matrix 'features' and a vector 'y' of target values
          row.names(features)=NULL

          # features with all values missing removed
          na_col = names ( which ( apply ( features, 2, function(x) all ( is.na ( x ) ) ) ) )
          features = features[,!names(features) %in% na_col]

          # features with infinite values removed
          inf_col = names ( which ( apply ( features, 2, function(x) any ( is.infinite ( x ) ) ) ) )
          features = features[,!names(features) %in% inf_col]

          # features with zero variance removed
          zero_var = names ( which ( apply ( features, 2, function(x) var(x, na.rm=T) ) == 0 ) )
          features = features[,!names(features) %in% zero_var]

          pp = NULL
          if (del_missing) {
            # needed if rows should be removed
            na_ids = apply ( features,1,function(x) any ( is.na ( x ) ) )
            features = features[!na_ids,]
            y = y[!na_ids]
            pp = preProcess(features, method=c("scale", "center"))
          } else {
            # Use imputation if NA's random (only then!)
            pp = preProcess(features, method=c("scale", "center", "knnImpute"))
          }
          features = predict(pp, features)

          # features with nan values removed (sometimes preProcess return NaN values)
          nan_col = names ( which ( apply ( features, 2, function(x) any ( is.nan ( x ) ) ) ) )
          features = features[,!names(features) %in% nan_col]

          # determine subsets
          subsets = dim(features)[2]*c(0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7)
          #subsets = dim(features)[2]*c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
          #subsets = c(2,3,4,5,7,10,subsets)
          #subsets = c(2,3,4,5,7,10,13,16,19,22,25,28,30)
          subsets = unique(sort(round(subsets))) 
          subsets = subsets[subsets<=dim(features)[2]]
          subsets = subsets[subsets>1] 

          # Recursive feature elimination
          rfProfile = rfe( x=features, y=y, rfeControl=rfeControl(functions=rfFuncs, number=150), sizes=subsets)

          # read existing dataset and select most useful features
          csv=feats[,c("SMILES", rfProfile$optVariables)]
          write.csv(x=csv,file=f_fds_r, row.names=F, quote=F, na='')
        EOR
        r_result_file
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
        params[:compound].lookup(params[:features], params[:feature_dataset_uri], params[:pc_type], params[:lib], params[:subjectid])
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

class Array
  # collect method extended for parallel processing.
  # Note: assign return value as: ans = arr.pcollect(n) { |obj| ... }
  # @param n the number of processes to spawn (default: unlimited)
  def pcollect(n = nil)
    nproc = 0
    result = collect do |*a|
      r, w = IO.pipe
      fork do
        r.close
        w.write( Marshal.dump( yield(*a) ) )
      end
      if n and (nproc+=1) >= n
        Process.wait ; nproc -= 1
      end
      [ w.close, r ].last
    end
    Process.waitall
    result.collect{|r| Marshal.load [ r.read, r.close ].first}
  end
end

