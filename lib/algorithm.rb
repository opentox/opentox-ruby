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
    module Fminer
      include Algorithm

      # Backbone Refinement Class mining (http://bbrc.maunz.de/)
      class BBRC
        include Fminer
        # Initialize bbrc algorithm
        def initialize(subjectid=nil)
          super File.join(CONFIG[:services]["opentox-algorithm"], "fminer/bbrc")
          load_metadata(subjectid)
        end
      end

      # LAtent STructure Pattern Mining (http://last-pm.maunz.de)
      class LAST
        include Fminer
        # Initialize last algorithm
        def initialize(subjectid=nil)
          super File.join(CONFIG[:services]["opentox-algorithm"], "fminer/last")
          load_metadata(subjectid)
        end
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

      # Classification with majority vote from neighbors weighted by similarity
      # @param [Array] neighbors, each neighbor is a hash with keys `:similarity, :activity`
      # @param [optional] params Ignored (only for compatibility with local_svm_regression)
      # @return [Hash] Hash with keys `:prediction, :confidence`
      def self.weighted_majority_vote(neighbors,params={}, props=nil)
        conf = 0.0
        confidence = 0.0
        neighbors.each do |neighbor|
          case neighbor[:activity].to_s
          when 'true'
            conf += Algorithm.gauss(neighbor[:similarity])
          when 'false'
            conf -= Algorithm.gauss(neighbor[:similarity])
          end
        end
        if conf > 0.0
          prediction = true
        elsif conf < 0.0
          prediction = false
        else
          prediction = nil
        end
        confidence = conf/neighbors.size if neighbors.size > 0
        {:prediction => prediction, :confidence => confidence.abs}
      end

      # Local support vector regression from neighbors 
      # @param [Array] neighbors, each neighbor is a hash with keys `:similarity, :activity, :features`
      # @param [Hash] params Keys `:similarity_algorithm,:p_values` are required
      # @return [Hash] Hash with keys `:prediction, :confidence`
      def self.local_svm_regression(neighbors, params)
        take_logs=true
        neighbors.each do |n| 
          if (! n[:activity].nil?) && (n[:activity].to_f < 0.0)
            take_logs = false
          end
        end
        acts = neighbors.collect do |n|
          act = n[:activity] 
          take_logs ? Math.log10(act.to_f) : act.to_f
        end # activities of neighbors for supervised learning

        sims = neighbors.collect{ |n| Algorithm.gauss(n[:similarity]) } # similarity values btwn q and nbors
        begin
          prediction = (props.nil? ? local_svm(neighbors, acts, sims, "nu-svr", params) : local_svm_prop(props, acts, "nu-svr", params))
          prediction = (take_logs ? 10**(prediction.to_f) : prediction.to_f)
          LOGGER.debug "Prediction is: '" + prediction.to_s + "'."
        rescue Exception => e
          LOGGER.debug "#{e.class}: #{e.message} #{e.backtrace}"
        end

        conf = sims.inject{|sum,x| sum + x }
        confidence = conf/neighbors.size if neighbors.size > 0
        {:prediction => prediction, :confidence => confidence}
        
      end

      # Local support vector classification from neighbors 
      # @param [Array] neighbors, each neighbor is a hash with keys `:similarity, :activity, :features`
      # @param [Hash] params Keys `:similarity_algorithm,:p_values` are required
      # @param [Array] props, propositionalization of neighbors and query structure e.g. [ Array_for_q, two-nested-Arrays_for_n ]
      # @return [Hash] Hash with keys `:prediction, :confidence`
      def self.local_svm_classification(neighbors, params, props=nil)
        acts = neighbors.collect do |n|
          act = n[:activity]
        end # activities of neighbors for supervised learning
        acts_f = acts.collect {|v| v == true ? 1.0 : 0.0}
        sims = neighbors.collect{ |n| Algorithm.gauss(n[:similarity]) } # similarity values btwn q and nbors
        begin 
          prediction = (props.nil? ? local_svm(neighbors, acts_f, sims, "C-bsvc", params) : local_svm_prop(props, acts_f, "C-bsvc", params))
          LOGGER.debug "Prediction is: '" + prediction.to_s + "'."
        rescue Exception => e
          LOGGER.debug "#{e.class}: #{e.message} #{e.backtrace}"
        end

        conf = sims.inject{|sum,x| sum + x }
        confidence = conf/neighbors.size if neighbors.size > 0
        {:prediction => prediction, :confidence => confidence}
        
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
                prediction = (@r.p.to_f == 1.0 ? true : false)
              end
              @r.quit # free R
            rescue Exception => e
              LOGGER.debug "#{e.class}: #{e.message} #{e.backtrace}"
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

          n_prop = props[0] # is a matrix, i.e. two nested Arrays.
          q_prop = props[1] # is an Array.

          #neighbor_matches = neighbors.collect{ |n| n[:features] } # URIs of matches
          #gram_matrix = [] # square matrix of similarities between neighbors; implements weighted tanimoto kernel
          if n_prop.size == 0
            raise "No neighbors found."
          else
            # gram matrix
            #(0..(neighbor_matches.length-1)).each do |i|
            #  gram_matrix[i] = [] unless gram_matrix[i]
            #  # upper triangle
            #  ((i+1)..(neighbor_matches.length-1)).each do |j|
            #    sim = eval("#{params[:similarity_algorithm]}(neighbor_matches[i], neighbor_matches[j], params[:p_values])")
            #    gram_matrix[i][j] = Algorithm.gauss(sim)
            #    gram_matrix[j] = [] unless gram_matrix[j] 
            #    gram_matrix[j][i] = gram_matrix[i][j] # lower triangle
            #  end
            #  gram_matrix[i][i] = 1.0
            #end

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
                prediction = (@r.p.to_f == 1.0 ? true : false)
              end
              @r.quit # free R
            rescue Exception => e
              LOGGER.debug "#{e.class}: #{e.message} #{e.backtrace}"
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

  end
end
