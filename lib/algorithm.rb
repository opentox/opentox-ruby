# R integration
# workaround to initialize R non-interactively (former rinruby versions did this by default)
# avoids compiling R with X
R = nil
require "rinruby" 

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
      
      def check_params(params,per_mil)
        raise OpenTox::NotFoundError.new "Please submit a dataset_uri." unless params[:dataset_uri] and  !params[:dataset_uri].nil?
        raise OpenTox::NotFoundError.new "Please submit a prediction_feature." unless params[:prediction_feature] and  !params[:prediction_feature].nil?
        @prediction_feature = OpenTox::Feature.find params[:prediction_feature], @subjectid
        @training_dataset = OpenTox::Dataset.find "#{params[:dataset_uri]}", @subjectid
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
      # @return [Float] (Weighted) tanimoto similarity
      def self.tanimoto(features_a,features_b,weights=nil)
        common_features = features_a & features_b
        all_features = (features_a + features_b).uniq
        common_p_sum = 0.0
        if common_features.size > 0
          if weights
            common_features.each{|f| common_p_sum += Algorithm.gauss(weights[f])}
            all_p_sum = 0.0
            all_features.each{|f| all_p_sum += Algorithm.gauss(weights[f])}
            common_p_sum/all_p_sum
          else
            common_features.to_f/all_features
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

    module Neighbors

      # Local multi-linear regression (MLR) prediction from neighbors. 
      # Uses propositionalized setting.
      # @param [Array] neighbors, each neighbor is a hash with keys `:similarity, :activity, :features`
      # @param [Hash] params Keys `:similarity_algorithm,:p_values` are required
      # @param [Array] props, propositionalization of neighbors and query structure e.g. [ Array_for_q, two-nested-Arrays_for_n ]
      # @return [Numeric] A prediction value.
      def self.local_mlr_prop(neighbors, params, props)

        raise "No neighbors found." unless neighbors.size>0
        begin

          acts = neighbors.collect do |n|
            act = n[:activity] 
            act.to_f
          end # activities of neighbors for supervised learning

          LOGGER.debug "Local MLR (Propositionalization / GSL)."
          n_prop = props[0] # is a matrix, i.e. two nested Arrays.
          q_prop = props[1] # is an Array.
          n_prop_x_size = n_prop[0].size
          n_prop_y_size = n_prop.size

          n_prop.flatten!
          y_x_rel = n_prop_y_size.to_f / n_prop_x_size
          repeat_factor = (1/y_x_rel).ceil
          n_prop_tmp = Array.new ; repeat_factor.times { n_prop_tmp.concat n_prop } ; n_prop = n_prop_tmp
          acts_tmp = Array.new ; repeat_factor.times { acts_tmp.concat acts } ; acts = acts_tmp

          LOGGER.debug "Setting GSL data ..."
          # set data
          prop_matrix = GSL::Matrix[n_prop, n_prop_y_size * repeat_factor, n_prop_x_size]
          y = GSL::Vector[acts]
          q_prop = GSL::Vector[q_prop]

          # model + support vectors
          LOGGER.debug "Creating MLR model ..."
          work = GSL::MultiFit::Workspace.alloc(n_prop_y_size * repeat_factor, n_prop_x_size)
          c, cov, chisq, status = GSL::MultiFit::linear(prop_matrix, y, work)
          LOGGER.debug "Predicting ..."
          prediction = GSL::MultiFit::linear_est(q_prop, c, cov)[0]
          LOGGER.debug "Prediction is: '" + prediction.to_s + "'."

          sims = neighbors.collect{ |n| Algorithm.gauss(n[:similarity]) } # similarity values btwn q and nbors
          conf = sims.inject{|sum,x| sum + x }
          confidence = conf/neighbors.size if neighbors.size > 0
          {:prediction => prediction, :confidence => confidence}

        rescue Exception => e
          LOGGER.debug "#{e.class}: #{e.message}"
          puts "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
        end


    end

      # Classification with majority vote from neighbors weighted by similarity
      # @param [Array] neighbors, each neighbor is a hash with keys `:similarity, :activity`
      # @param [optional] params Ignored (only for compatibility with local_svm_regression)
      # @return [Hash] Hash with keys `:prediction, :confidence`
      def self.weighted_majority_vote(neighbors,params={}, props=nil)
        neighbor_contribution = 0.0
        confidence_sum = 0.0
        confidence = 0.0
        prediction = nil
        positive_map_value= nil
        negative_map_value= nil

        neighbors.each do |neighbor|
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
            prediction = 2 unless neighbors.size==0
          elsif confidence_sum < 0.0
            prediction = 1 unless neighbors.size==0
          end
        else 
          prediction = (neighbor_contribution/confidence_sum).round  unless neighbors.size==0  # AM: new multinomial prediction
        end 

        confidence = confidence_sum/neighbors.size if neighbors.size > 0
        return {:prediction => prediction, :confidence => confidence.abs}
      end

      # Local support vector regression from neighbors 
      # @param [Array] neighbors, each neighbor is a hash with keys `:similarity, :activity, :features`
      # @param [Hash] params Keys `:similarity_algorithm,:p_values` are required
      # @return [Hash] Hash with keys `:prediction, :confidence`
      def self.local_svm_regression(neighbors, params, props=nil)

        raise "No neighbors found." unless neighbors.size>0
        begin
          acts = neighbors.collect{ |n| n[:activity].to_f }
          sims = neighbors.collect{ |n| Algorithm.gauss(n[:similarity]) }
          offset = 1.0 - acts.minmax[0] # offset to min element
          offset = -1.0 * offset if offset>0.0 
          inverter = OpenTox::Algorithm::Transform::Inverter.new(acts)
          prediction = (props.nil? ? local_svm(neighbors, inverter.values, sims, "nu-svr", params) : local_svm_prop(props, inverter.values, "nu-svr", params))
          prediction = inverter.back_transform([prediction])[0]
          LOGGER.debug "Prediction is: '" + prediction.to_s + "'."
          conf = sims.inject{|sum,x| sum + x }
          confidence = conf/neighbors.size
          {:prediction => prediction, :confidence => confidence}
        rescue Exception => e
          LOGGER.debug "#{e.class}: #{e.message}"
          puts "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
        end
        
      end

      # Local support vector classification from neighbors 
      # @param [Array] neighbors, each neighbor is a hash with keys `:similarity, :activity, :features`
      # @param [Hash] params Keys `:similarity_algorithm,:p_values` are required
      # @param [Array] props, propositionalization of neighbors and query structure e.g. [ Array_for_q, two-nested-Arrays_for_n ]
      # @return [Hash] Hash with keys `:prediction, :confidence`
      def self.local_svm_classification(neighbors, params, props=nil)

        raise "No neighbors found." unless neighbors.size>0
        begin 
          acts = neighbors.collect { |n| act = n[:activity] }
          acts_f = acts
          sims = neighbors.collect{ |n| Algorithm.gauss(n[:similarity]) } # similarity values btwn q and nbors
          prediction = (props.nil? ? local_svm(neighbors, acts_f, sims, "C-bsvc", params) : local_svm_prop(props, acts_f, "C-bsvc", params))
          LOGGER.debug "Prediction is: '" + prediction.to_s + "'."
          conf = sims.inject{|sum,x| sum + x }
          confidence = conf/neighbors.size if neighbors.size > 0
          {:prediction => prediction, :confidence => confidence}
        rescue Exception => e
          LOGGER.debug "#{e.class}: #{e.message}"
          puts "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
        end
        
      end


      # Local support vector prediction from neighbors. 
      # Uses pre-defined Kernel Matrix.
      # Not to be called directly (use local_svm_regression or local_svm_classification).
      # @param [Array] neighbors, each neighbor is a hash with keys `:similarity, :activity, :features`
      # @param [Array] acts, activities for neighbors.
      # @param [Array] sims, similarities for neighbors.
      # @param [String] type, one of "nu-svr" (regression) or "C-bsvc" (classification).
      # @param [Hash] params Keys `:similarity_algorithm,:p_values` are required
      # @param [Array] props, propositionalization of neighbors and query structure e.g. [ Array_for_q, two-nested-Arrays_for_n ]
      # @return [Numeric] A prediction value.
      def self.local_svm(neighbors, acts, sims, type, params)
        LOGGER.debug "Local SVM (Weighted Tanimoto Kernel)."
        neighbor_matches = neighbors.collect{ |n| n[:features] } # URIs of matches
        gram_matrix = [] # square matrix of similarities between neighbors; implements weighted tanimoto kernel
        if neighbor_matches.size == 0
          raise "No neighbors found."
        else
          # gram matrix
          (0..(neighbor_matches.length-1)).each do |i|
            gram_matrix[i] = [] unless gram_matrix[i]
            # upper triangle
            ((i+1)..(neighbor_matches.length-1)).each do |j|
              sim = eval("#{params[:similarity_algorithm]}(neighbor_matches[i], neighbor_matches[j], params[:p_values])")
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
            puts "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
          end

        end
        prediction
      end

      # Local support vector prediction from neighbors. 
      # Uses propositionalized setting.
      # Not to be called directly (use local_svm_regression or local_svm_classification).
      # @param [Array] neighbors, each neighbor is a hash with keys `:similarity, :activity, :features`
      # @param [Array] acts, activities for neighbors.
      # @param [Array] props, propositionalization of neighbors and query structure e.g. [ Array_for_q, two-nested-Arrays_for_n ]
      # @param [String] type, one of "nu-svr" (regression) or "C-bsvc" (classification).
      # @param [Hash] params Keys `:similarity_algorithm,:p_values` are required
      # @return [Numeric] A prediction value.
      def self.local_svm_prop(props, acts, type, params)

          LOGGER.debug "Local SVM (Propositionalization / Kernlab Kernel)."
          n_prop = props[0] # is a matrix, i.e. two nested Arrays.
          q_prop = props[1] # is an Array.

          #neighbor_matches = neighbors.collect{ |n| n[:features] } # URIs of matches
          #gram_matrix = [] # square matrix of similarities between neighbors; implements weighted tanimoto kernel
          if n_prop.size == 0
            raise "No neighbors found."
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
              puts "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
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

    module Transform
      include Algorithm

      class Inverter # to improve normality conditions on a vector
        attr_accessor :values

        def initialize(values)
          @values=values
          raise "Cannot transform, values empty." if @values.size==0
          @values = @values.collect { |v| -1.0 * v }   # slide >1
          @offset = 1.0 - @values.minmax[0] 
          @offset = -1.0 * @offset if @offset>0.0 
          @values = @values.collect { |v| v - @offset }   # slide >1
          @values = @values.collect { |v| 1 / v }         # invert using sigmoidal function
        end

        def back_transform(values)
          values = values.collect { |v| 1 / v }
          values = values.collect { |v| v + @offset }
          values = values.collect { |v| -1.0 * v }
        end

      end
    end
    
    # Gauss kernel
    # @return [Float] 
    def self.gauss(x, sigma = 0.3) 
      d = 1.0 - x.to_f
      Math.exp(-(d*d)/(2*sigma*sigma))
    end
    
    # Median of an array
    # @param [Array] Array with values
    # @return [Float] Median
    def self.median(array)
      return nil if array.empty?
      array.sort!
      m_pos = array.size / 2
      return array.size % 2 == 1 ? array[m_pos] : (array[m_pos-1] + array[m_pos])/2
    end

    # Sum of an array for Numeric values
    # @param [Array] Array with values
    # @return [Integer] Sum of values
    def self.sum(array)
      array.inject{|s,x| s + x }
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
      minfreq = per_mil*training_dataset.compounds.size/1000 # AM sugg. 8-10 per mil for BBRC, 50 per mil for LAST
      minfreq = 2 unless minfreq > 2
      minfreq
    end

    # Effect calculation for classification
    # @param [Array] Array of occurrences per class in the form of Enumerables.
    # @param [Array] Array of database instance counts per class.
    def self.effect(occurrences, db_instances)
      max=0
      max_value=0
      nr_o = self.sum_size(occurrences)
      nr_db = self.sum(db_instances)

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


  end
end


