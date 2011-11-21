module OpenTox
  module Algorithm
   
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
      return (array.to_scale.variance_population == 0.0)
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

    module Neighbors

      # Calculate the propositionalization matrix (aka instantiation matrix) via fingerprints.
      # Same for the vector describing the query compound.
      # @param[Hash] Required keys: :neighbors, compound, features, nr_hits, fingerprints, p_values
      def self.get_props_fingerprints (params)
        matrix = Array.new
        begin 
          
          # neighbors
          params[:neighbors].each do |n|
            n = n[:compound]
            row = []
            row_good = true
            params[:features].each do |f|
              #if (!params[:fingerprints][n].nil?) && (params[:fingerprints][n].include?(f))
              #  row << params[:p_values][f] * params[:fingerprints][n][f]
              #else
              #  LOGGER.debug "Warning: Neighbor with missing values skipped." if row_good
              #  row_good = false
              #end
              if ! params[:fingerprints][n].nil? 
                row << (params[:fingerprints][n].include?(f) ? (params[:p_values][f] * params[:fingerprints][n][f]) : 0.0)
              else
                row << 0.0
              end
            end
            matrix << row if row_good
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
          LOGGER.debug "get_props_fingerprints failed with '" + $! + "'"
        end
        [ matrix, row ]
      end
      
      
      # Get X and Y size of a nested Array (Matrix)
      # @param [Array] Two-dimensional ruby array (matrix) with X and Y size > 0
      # @return [Arrray] X and Y size of the matrix
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

    end

  end

end

