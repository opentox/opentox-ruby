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
        return @feature_type if @feature_type

        # dynamically perform restcalls if necessary
        load_metadata(subjectid) if @metadata==nil or @metadata.size==0 or (@metadata.size==1 && @metadata.values[0]==@uri)
        algorithm = OpenTox::Algorithm::Generic.find(@metadata[OT.algorithm], subjectid)
        algorithm_title = algorithm ? algorithm.metadata[DC.title] : nil
        algorithm_type = algorithm ? algorithm.metadata[RDF.type] : nil
        dependent_variable = OpenTox::Feature.find( @metadata[OT.dependentVariables],subjectid )
        dependent_variable_type = dependent_variable ? dependent_variable.feature_type : nil
        type_indicators = [dependent_variable_type, @metadata[RDF.type], @metadata[DC.title], @uri, algorithm_type, algorithm_title].flatten
        type_indicators.each do |type|
          case type
          when /(?i)classification/
            @feature_type =  "classification"
            break
          when /(?i)regression/
            @feature_type = "regression"
          end
        end
        raise "unknown model "+type_indicators.inspect unless @feature_type
        @feature_type
      end

    end

    # Lazy Structure Activity Relationship class
    class Lazar

      include Model
      include Algorithm

      attr_accessor :compound, :prediction_dataset, :features, :effects, :activities, :p_values, :fingerprints, :feature_calculation_algorithm, :similarity_algorithm, :prediction_algorithm, :min_sim, :subjectid

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

        @feature_calculation_algorithm = "Substructure.match"
        @similarity_algorithm = "Similarity.tanimoto"
        @prediction_algorithm = "Neighbors.weighted_majority_vote"

        @min_sim = 0.3

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
      def self.create(params)
        subjectid = params[:subjectid]
        lazar_algorithm = OpenTox::Algorithm::Generic.new File.join( CONFIG[:services]["opentox-algorithm"],"lazar")
        model_uri = lazar_algorithm.run(params)
        OpenTox::Model::Lazar.find(model_uri, subjectid)      
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
        @prediction_dataset.save(subjectid)
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
            # TODO: fix dependentVariable
            DC.title => URI.decode(File.basename( @metadata[OT.dependentVariables] )),
            OT.parameters => [{DC.title => "compound_uri", OT.paramValue => compound_uri}]
          } )
        end

        return @prediction_dataset if database_activity(subjectid)


       if metadata[RDF.type] == [OTA.ClassificationLazySingleTarget]
          # AM: Balancing, see http://www.maunz.de/wordpress/opentox/2011/balanced-lazar
          l = Array.new # larger 
          s = Array.new # smaller fraction
          @fingerprints.each do |training_compound,training_features|
            @activities[training_compound].each do |act|
              case act.to_s
              when "false" 
                l << training_compound
              when "true"  
                s << training_compound
              else
                LOGGER.warn "BLAZAR: Activity #{act.to_s} should not be reached."
              end
            end
          end
          if s.size > l.size then 
            l,s = s,l # happy swapping
            LOGGER.info "BLAZAR: |s|=#{s.size}, |l|=#{l.size}."
          end
          # determine ratio
          modulo = l.size.divmod(s.size)# modulo[0]=ratio, modulo[1]=rest
          LOGGER.info "BLAZAR: Balance: #{modulo[0]}, rest #{modulo[1]}."

          # AM: Balanced predictions
          addon = (modulo[1].to_f/modulo[0]).ceil # what will be added in each round 
          slack = modulo[1].divmod(addon)[1] # what remains for the last round
          position = 0
          predictions = Array.new

          prediction_best=nil
          neighbors_best=nil

          begin
          for i in 1..modulo[0] do
            (i == modulo[0]) && (slack>0) ? lr_size = s.size + slack : lr_size = s.size + addon  # determine fraction
            LOGGER.info "BLAZAR: Neighbors round #{i}: #{position} + #{lr_size}."
            neighbors_balanced(s, l, position, lr_size) # get ratio fraction of larger part
            prediction = eval("#{@prediction_algorithm}(@neighbors,{:similarity_algorithm => @similarity_algorithm, :p_values => @p_values})")
            if prediction_best.nil? || prediction[:confidence].abs > prediction_best[:confidence].abs 
              prediction_best=prediction 
              neighbors_best=@neighbors
            end
            position = position + lr_size
          end
          rescue Exception => e
            LOGGER.error "BLAZAR failed in prediction: "+e.class.to_s+": "+e.message
          end

          prediction=prediction_best
          @neighbors=neighbors_best
          ### END AM balanced predictions

        else # regression case: no balancing
          neighbors
          prediction = eval("#{@prediction_algorithm}(@neighbors,{:similarity_algorithm => @similarity_algorithm, :p_values => @p_values})")
        end
        
       # TODO: reasonable feature name
        #prediction_feature_uri = File.join( @prediction_dataset.uri, "feature", "prediction", File.basename(@metadata[OT.dependentVariables]),@prediction_dataset.compounds.size.to_s)
        value_feature_uri = File.join( @prediction_dataset.uri, "feature", "prediction", File.basename(@metadata[OT.dependentVariables]),"value")
        confidence_feature_uri = File.join( @prediction_dataset.uri, "feature", "prediction", File.basename(@metadata[OT.dependentVariables]),"confidence")

        prediction_feature_uris = {value_feature_uri => prediction[:prediction], confidence_feature_uri => prediction[:confidence]}
        #prediction_feature_uris[value_feature_uri] = "No similar compounds in training dataset." if @neighbors.size == 0 or prediction[:prediction].nil?
        prediction_feature_uris[value_feature_uri] = nil if @neighbors.size == 0 or prediction[:prediction].nil?
        
        #@prediction_dataset.metadata[OT.dependentVariables] = prediction_feature_uri
        @prediction_dataset.metadata[OT.dependentVariables] = @metadata[OT.dependentVariables]

=begin
        if @neighbors.size == 0
          prediction_feature_uris.each do |prediction_feature_uri,value|
            @prediction_dataset.add_feature(prediction_feature_uri, {
              RDF.type => [OT.MeasuredFeature],
              OT.hasSource => @uri,
              DC.creator => @uri,
              DC.title => URI.decode(File.basename( @metadata[OT.dependentVariables] )),
              OT.error => "No similar compounds in training dataset.",
              #OT.parameters => [{DC.title => "compound_uri", OT.paramValue => compound_uri}]
            })
            @prediction_dataset.add @compound.uri, prediction_feature_uri, value
          end

        else
=end
        prediction_feature_uris.each do |prediction_feature_uri,value|
          @prediction_dataset.metadata[OT.predictedVariables] = [] unless @prediction_dataset.metadata[OT.predictedVariables] 
          @prediction_dataset.metadata[OT.predictedVariables] << prediction_feature_uri
          @prediction_dataset.add_feature(prediction_feature_uri, {
            RDF.type => [OT.ModelPrediction],
            OT.hasSource => @uri,
            DC.creator => @uri,
            DC.title => URI.decode(File.basename( @metadata[OT.dependentVariables] )),
            # TODO: factor information to value
          })
            #OT.prediction => prediction[:prediction],
            #OT.confidence => prediction[:confidence],
            #OT.parameters => [{DC.title => "compound_uri", OT.paramValue => compound_uri}]
            @prediction_dataset.add @compound.uri, prediction_feature_uri, value
        end

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
        #end

        @prediction_dataset.save(subjectid)
        @prediction_dataset
      end

      # Calculate the propositionalization matrix aka instantiation matrix (0/1 entries for features)
#      def get_prop_matrix
#        matrix = Array.new
#        begin 
#        @neighbors.each do |n|
#          row = []
#          @features.each do |f|
#            row << @fingerprints[n].include?(f) ? 0.0 : @p_values[f]
#          end
#          matrix << row
#        end
#        rescue Exception => e
#          LOGGER.debug "get_prop_matrix failed with '" + $! + "'"
#        end
#        matrix
#      end

      # Find neighbors and store them as object variable, access only a subset of compounds for that.
      def neighbors_balanced(s, l, start, offset)
        @compound_features = eval("#{@feature_calculation_algorithm}(@compound,@features)") if @feature_calculation_algorithm
        @neighbors = []
          [ l[start, offset ] , s ].flatten.each do |training_compound| # AM: access only a balanced subset
            training_features = @fingerprints[training_compound]
            add_neighbor training_features, training_compound
        end

      end

      # Find neighbors and store them as object variable, access all compounds for that.
      def neighbors
        @compound_features = eval("#{@feature_calculation_algorithm}(@compound,@features)") if @feature_calculation_algorithm
        @neighbors = []
        @fingerprints.each do |training_compound,training_features| # AM: access all compounds
           add_neighbor training_features, training_compound
        end
      end

      # Adds a neighbor to @neighbors if it passes the similarity threshold.
      def add_neighbor(training_features, training_compound)
        sim = eval("#{@similarity_algorithm}(@compound_features,training_features,@p_values)")
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
          @activities[@compound.uri].each { |act| @prediction_dataset.add @compound.uri, @metadata[OT.dependentVariables], act }
          @prediction_dataset.add_metadata(OT.hasSource => @metadata[OT.trainingDataset])
          @prediction_dataset.save(subjectid)
          true
        else
          false
        end
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
