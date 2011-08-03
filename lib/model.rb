module OpenTox

  module Model

    include OpenTox

    # Run a model with parameters
    # @param [Hash] params Parameters for OpenTox model
    # @param [optional,OpenTox::Task] waiting_task (can be a OpenTox::Subtask as well), progress is updated accordingly
    # @return [text/uri-list] Task or resource URI
    def run( params, accept_header=nil, waiting_task=nil )
      unless accept_header
        if CONFIG[:yaml_hosts].include?(URI.parse(@uri).host)
          accept_header = 'application/x-yaml' 
        else
          accept_header = 'application/rdf+xml'
        end
      end
      LOGGER.info "running model "+@uri.to_s+", params: "+params.inspect+", accept: "+accept_header.to_s
      RestClientWrapper.post(@uri,params,{:accept => accept_header},waiting_task).to_s
    end

    # Generic OpenTox model class for all API compliant services
    class Generic
      include Model

      # Find Generic Opentox Model via URI, and loads metadata, could raise NotFound/NotAuthorized error 
      # @param [String] uri Model URI
      # @return [OpenTox::Model::Generic] Model instance
      def self.find(uri,subjectid=nil)
        return nil unless uri
        model = Generic.new(uri)
        model.load_metadata(subjectid)
        raise "could not load model metadata '"+uri.to_s+"'" if model.metadata==nil or model.metadata.size==0
        model
      end

      # provides feature type, possible types are "regression" or "classification"
      # @return [String] feature type, "unknown" if type could not be estimated
      def feature_type(subjectid=nil)
        unless @feature_type
          load_predicted_variables( subjectid ) unless @predicted_variable
          @feature_type = OpenTox::Feature.find( @predicted_variable, subjectid ).feature_type
        end
        @feature_type
      end
    
      def predicted_variable( subjectid )
        load_predicted_variables( subjectid ) unless @predicted_variable
        @predicted_variable
      end

      def predicted_variables( subjectid )
        load_predicted_variables( subjectid, false ) unless @predicted_variables
        @predicted_variables
      end

      def predicted_confidence( subjectid )
        load_predicted_variables( subjectid ) unless @predicted_confidence
        @predicted_confidence
      end
  
      private
      def load_predicted_variables( subjectid=nil, use_confidence=true )
        load_metadata(subjectid) if @metadata==nil or @metadata.size==0 or (@metadata.size==1 && @metadata.values[0]==@uri)
        if @metadata[OT.predictedVariables]
          predictedVariables = @metadata[OT.predictedVariables]
          if predictedVariables.is_a?(Array)
            if (predictedVariables.size==1)
              @predicted_variable = predictedVariables[0]
            elsif (predictedVariables.size>=2)
              # PENDING identify confidence
              if use_confidence
                conf_index = -1
                predictedVariables.size.times do |i|
                  f = OpenTox::Feature.find(predictedVariables[i], subjectid)
                  conf_index = i if f.metadata[DC.title]=~/(?i)confidence/
                end
                raise "could not estimate predicted variable from model: '"+uri.to_s+
                  "', number of predicted-variables==2, but no confidence found" if conf_index==-1
              end
              if (predictedVariables.size==2) && use_confidence
                @predicted_variable = predictedVariables[1-conf_index]
                @predicted_confidence = predictedVariables[conf_index]
              else
                @predicted_variables = predictedVariables
              end
            else
              raise "could not estimate predicted variable from model: '"+uri.to_s+"', number of predicted-variables == 0"  
            end
          else
            raise "could not estimate predicted variable from model: '"+uri.to_s+"', predicted-variables is no array"
          end        
        end
        raise "could not estimate predicted variable from model: '"+uri.to_s+"'" unless (@predicted_variable || @predicted_variables)
      end
    end

    # Lazy Structure Activity Relationship class
    class Lazar

      include Algorithm
      include Model

      attr_accessor :compound, :prediction_dataset, :features, :effects, :activities, :p_values, :fingerprints, :feature_calculation_algorithm, :similarity_algorithm, :prediction_algorithm, :min_sim, :subjectid, :prop_kernel, :value_map, :nr_hits, :transform, :conf_stdev

      def initialize(uri=nil)

        if uri
          super uri
        else
          super CONFIG[:services]["opentox-model"]
        end

        @metadata[OT.algorithm] = File.join(CONFIG[:services]["opentox-algorithm"],"lazar")

        @features = []
        @effects = {}
        @activities = {}
        @p_values = {}
        @fingerprints = {}
        @value_map = {}

        @feature_calculation_algorithm = "Substructure.match"
        @similarity_algorithm = "Similarity.tanimoto"
        @prediction_algorithm = "Neighbors.weighted_majority_vote"
        
        @nr_hits = false
        @min_sim = 0.3
        @prop_kernel = false
        @transform = { "class" => "NOP"  }
        @conf_stdev = false

      end

      # Get URIs of all lazar models
      # @return [Array] List of lazar model URIs
      def self.all(subjectid=nil)
        RestClientWrapper.get(CONFIG[:services]["opentox-model"], :subjectid => subjectid).to_s.split("\n")
      end

      # Find a lazar model
      # @param [String] uri Model URI
      # @return [OpenTox::Model::Lazar] lazar model
      def self.find(uri, subjectid=nil)
        YAML.load RestClientWrapper.get(uri,{:accept => 'application/x-yaml', :subjectid => subjectid})
      end

      # Create a new lazar model
      # @param [optional,Hash] params Parameters for the lazar algorithm (OpenTox::Algorithm::Lazar)
      # @return [OpenTox::Model::Lazar] lazar model
      def self.create(params, waiting_task=nil )
        subjectid = params[:subjectid]
        lazar_algorithm = OpenTox::Algorithm::Generic.new File.join( CONFIG[:services]["opentox-algorithm"],"lazar")
        model_uri = lazar_algorithm.run(params, waiting_task)
        OpenTox::Model::Lazar.find(model_uri, subjectid)      
      end

      def run( params, accept_header=nil, waiting_task=nil )
      unless accept_header
        if CONFIG[:yaml_hosts].include?(URI.parse(@uri).host)
          accept_header = 'application/x-yaml' 
        else
          accept_header = 'application/rdf+xml'
        end
      end
      LOGGER.info "running model "+@uri.to_s+", params: "+params.inspect+", accept: "+accept_header.to_s
      RestClientWrapper.post(@uri,params,{:accept => accept_header},waiting_task).to_s
      end

      # Get a parameter value
      # @param [String] param Parameter name
      # @return [String] Parameter value
      def parameter(param)
        @metadata[OT.parameters].collect{|p| p[OT.paramValue] if p[DC.title] == param}.compact.first
      end

      # Predict a dataset
      # @param [String] dataset_uri Dataset URI
      # @param [optional,subjectid] 
      # @param [optional,OpenTox::Task] waiting_task (can be a OpenTox::Subtask as well), progress is updated accordingly
      # @return [OpenTox::Dataset] Dataset with predictions
      def predict_dataset(dataset_uri, subjectid=nil, waiting_task=nil)
      
        @prediction_dataset = Dataset.create(CONFIG[:services]["opentox-dataset"], subjectid)
        @prediction_dataset.add_metadata({
          OT.hasSource => @uri,
          DC.creator => @uri,
          DC.title => URI.decode(File.basename( @metadata[OT.dependentVariables] )),
          OT.parameters => [{DC.title => "dataset_uri", OT.paramValue => dataset_uri}]
        })
        d = Dataset.new(dataset_uri,subjectid)
        d.load_compounds(subjectid)
        count = 0
        d.compounds.each do |compound_uri|
          begin
            predict(compound_uri,false,subjectid)
            count += 1
            waiting_task.progress( count/d.compounds.size.to_f*100.0 ) if waiting_task
          rescue => ex
            LOGGER.warn "prediction for compound "+compound_uri.to_s+" failed: "+ex.message
          end
        end
        #@prediction_dataset.save(subjectid)
        @prediction_dataset
      end

      # Predict a compound
      # @param [String] compound_uri Compound URI
      # @param [optinal,Boolean] verbose Verbose prediction (output includes neighbors and features)
      # @return [OpenTox::Dataset] Dataset with prediction
      def predict(compound_uri,verbose=false,subjectid=nil)

        @compound = Compound.new compound_uri
        features = {}

        unless @prediction_dataset
          @prediction_dataset = Dataset.create(CONFIG[:services]["opentox-dataset"], subjectid)
          @prediction_dataset.add_metadata( {
            OT.hasSource => @uri,
            DC.creator => @uri,
            DC.title => URI.decode(File.basename( @metadata[OT.dependentVariables] )),
            OT.parameters => [{DC.title => "compound_uri", OT.paramValue => compound_uri}]
          } )
        end

        unless database_activity(subjectid) # adds database activity to @prediction_dataset

          neighbors
          prediction = eval("#{@prediction_algorithm} ( { :neighbors => @neighbors, 
                                                          :compound => @compound,
                                                          :features => @features, 
                                                          :p_values => @p_values, 
                                                          :fingerprints => @fingerprints,
                                                          :similarity_algorithm => @similarity_algorithm, 
                                                          :prop_kernel => @prop_kernel,
                                                          :value_map => @value_map,
                                                          :nr_hits => @nr_hits,
                                                          :conf_stdev => @conf_stdev,
                                                          :transform => @transform } ) ")

          value_feature_uri = File.join( @uri, "predicted", "value")
          confidence_feature_uri = File.join( @uri, "predicted", "confidence")

          @prediction_dataset.metadata[OT.dependentVariables] = @metadata[OT.dependentVariables] unless @prediction_dataset.metadata[OT.dependentVariables] 
          @prediction_dataset.metadata[OT.predictedVariables] = [value_feature_uri, confidence_feature_uri] unless @prediction_dataset.metadata[OT.predictedVariables] 

          if OpenTox::Feature.find(metadata[OT.dependentVariables], subjectid).feature_type == "classification"
            @prediction_dataset.add @compound.uri, value_feature_uri, @value_map[prediction[:prediction]]
          else
            @prediction_dataset.add @compound.uri, value_feature_uri, prediction[:prediction]
          end
          @prediction_dataset.add @compound.uri, confidence_feature_uri, prediction[:confidence]
          @prediction_dataset.features[value_feature_uri][DC.title] = @prediction_dataset.metadata[DC.title]
          @prediction_dataset.features[confidence_feature_uri][DC.title] = "Confidence"

          if verbose
            if @feature_calculation_algorithm == "Substructure.match"
              f = 0
              @compound_features.each do |feature|
                feature_uri = File.join( @prediction_dataset.uri, "feature", "descriptor", f.to_s)
                features[feature] = feature_uri
                @prediction_dataset.add_feature(feature_uri, {
                  RDF.type => [OT.Substructure],
                  OT.smarts => feature,
                  OT.pValue => @p_values[feature],
                  OT.effect => @effects[feature]
                })
                @prediction_dataset.add @compound.uri, feature_uri, true
                f+=1
              end
            else
              @compound_features.each do |feature|
                features[feature] = feature
                @prediction_dataset.add @compound.uri, feature, true
              end
            end
            n = 0
            @neighbors.each do |neighbor|
              neighbor_uri = File.join( @prediction_dataset.uri, "feature", "neighbor", n.to_s )
              @prediction_dataset.add_feature(neighbor_uri, {
                OT.compound => neighbor[:compound],
                OT.similarity => neighbor[:similarity],
                OT.measuredActivity => neighbor[:activity],
                RDF.type => [OT.Neighbor]
              })
              @prediction_dataset.add @compound.uri, neighbor_uri, true
              f = 0 unless f
              neighbor[:features].each do |feature|
                if @feature_calculation_algorithm == "Substructure.match"
                  feature_uri = File.join( @prediction_dataset.uri, "feature", "descriptor", f.to_s) unless feature_uri = features[feature]
                else
                  feature_uri = feature
                end
                @prediction_dataset.add neighbor[:compound], feature_uri, true
                unless features.has_key? feature
                  features[feature] = feature_uri
                  @prediction_dataset.add_feature(feature_uri, {
                    RDF.type => [OT.Substructure],
                    OT.smarts => feature,
                    OT.pValue => @p_values[feature],
                    OT.effect => @effects[feature]
                  })
                  f+=1
                end
              end
              n+=1
            end
          end
        end

        @prediction_dataset.save(subjectid)
        @prediction_dataset
      end

      

      # Find neighbors and store them as object variable, access all compounds for that.
      def neighbors
        @compound_features = eval("#{@feature_calculation_algorithm}(@compound,@features)") if @feature_calculation_algorithm
        @neighbors = []
        @fingerprints.keys.each do |training_compound| # AM: access all compounds
          add_neighbor @fingerprints[training_compound].keys, training_compound
        end
      end

      # Adds a neighbor to @neighbors if it passes the similarity threshold.
      def add_neighbor(training_features, training_compound)
        compound_features_hits = {}
        training_compound_features_hits = {}
        if @nr_hits
          compound_features_hits = @compound.match_hits(@compound_features)
          training_compound_features_hits = @fingerprints[training_compound]
          #LOGGER.debug "dv ------------ training_compound_features_hits:#{training_compound_features_hits.class}  #{training_compound_features_hits}"
        end
        params = {}
        params[:nr_hits] = @nr_hits
        params[:compound_features_hits] = compound_features_hits
        params[:training_compound_features_hits] = training_compound_features_hits

        sim = eval("#{@similarity_algorithm}(training_features, @compound_features, @p_values, params)")
        if sim > @min_sim
          @activities[training_compound].each do |act|
            @neighbors << {
              :compound => training_compound,
              :similarity => sim,
              :features => training_features,
              :activity => act
            }
          end
        end
      end

      # Find database activities and store them in @prediction_dataset
      # @return [Boolean] true if compound has databasse activities, false if not
      def database_activity(subjectid)
        if @activities[@compound.uri]
          @activities[@compound.uri].each { |act| @prediction_dataset.add @compound.uri, @metadata[OT.dependentVariables], @value_map[act] }
          @prediction_dataset.add_metadata(OT.hasSource => @metadata[OT.trainingDataset])
          @prediction_dataset.save(subjectid)
          true
        else
          false
        end
      end

      def prediction_features
        [prediction_value_feature,prediction_confidence_feature]
      end

      def prediction_value_feature
        dependent_uri = @metadata[OT.dependentVariables].first
        feature = OpenTox::Feature.new File.join( @uri, "predicted", "value")
        feature.add_metadata( {
          RDF.type => [OT.ModelPrediction],
          OT.hasSource => @uri,
          DC.creator => @uri,
          DC.title => URI.decode(File.basename( dependent_uri )),
          OWL.sameAs => dependent_uri
        })
        feature
      end

      def prediction_confidence_feature
        dependent_uri = @metadata[OT.dependentVariables].first
        feature = OpenTox::Feature.new File.join( @uri, "predicted", "confidence")
        feature.add_metadata( {
          RDF.type => [OT.ModelPrediction],
          OT.hasSource => @uri,
          DC.creator => @uri,
          DC.title => "#{URI.decode(File.basename( dependent_uri ))} confidence"
        })
        feature
      end

      # Save model at model service
      def save(subjectid)
        self.uri = RestClientWrapper.post(@uri,self.to_yaml,{:content_type =>  "application/x-yaml", :subjectid => subjectid})
      end

      # Delete model at model service
      def delete(subjectid)
        RestClientWrapper.delete(@uri, :subjectid => subjectid) unless @uri == CONFIG[:services]["opentox-model"]
      end

    end
  end
end
