require 'csv'


module OpenTox

  module Algorithm

    # Calculate physico-chemical descriptors.
    # @param[Hash] Required keys: :dataset_uri, :pc_type
    # @return[String] dataset uri

    def self.pc_descriptors(params)

      begin
        ds = OpenTox::Dataset.find(params[:dataset_uri])
        compounds = ds.compounds.collect
        ambit_result_uri = get_pc_descriptors( { :compounds => compounds, :pc_type => params[:pc_type] } )
        #ambit_result_uri = "http://apps.ideaconsult.net:8080/ambit2/dataset/987103?feature_uris[]=http%3A%2F%2Fapps.ideaconsult.net%3A8080%2Fambit2%2Ffeature%2F4276789&feature_uris[]=http%3A%2F%2Fapps.ideaconsult.net%3A8080%2Fambit2%2Fmodel%2F16%2Fpredicted" # for testing
        LOGGER.debug "Ambit result uri: '#{ambit_result_uri}'"
        csv_data = CSV.parse( OpenTox::RestClientWrapper.get(ambit_result_uri, {:accept => "text/csv"}) )

        index_ambit_uri = csv_data[0].index("Compound")
        csv_data.map {|i| i.delete_at(index_ambit_uri)}
        csv_data[0].each {|cell| cell.chomp!(" ")}
        
        parser = OpenTox::Parser::Spreadsheets.new
        ds = OpenTox::Dataset.new
        ds.save
        parser.dataset = ds
        ds = parser.load_csv(csv_data.collect{|r| r.join(",")}.join("\n"))
        ds.save
      rescue Exception => e
        LOGGER.debug "#{e.class}: #{e.message}"
        LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
      end

    end

    # Calculates PC descriptors via Ambit -- DO NOT OVERLOAD Ambit.
    # @param[Hash] Required keys: :compounds, :pc_type
    # @return[String] Ambit result uri
    def self.get_pc_descriptors(params)

      begin

        ambit_ds_service_uri = "http://apps.ideaconsult.net:8080/ambit2/dataset/"
        ambit_mopac_model_uri = "http://apps.ideaconsult.net:8080/ambit2/model/69632"
        descs = YAML::load_file( File.join(ENV['HOME'], ".opentox", "config", "ambit_descriptors.yaml") )
        descs_uris = []
        params[:pc_type] = "electronic,cpsa" if params[:pc_type].nil? # rescue missing pc_type
        types = params[:pc_type].split(",")
        descs.each { |uri, cat_name| 
          if types.include? cat_name[:category]
            descs_uris << uri
          end
        }
        if descs_uris.size == 0
          raise "Error! Empty set of descriptors. Did you supply one of [geometrical, topological, electronic, constitutional, hybrid, cpsa] ?"
        end
        #LOGGER.debug "Ambit descriptor URIs: #{descs_uris.join(", ")}"

        begin
          # Create SMI
          smiles = []
          params[:compounds].each do |n|
            smiles << OpenTox::Compound.new(n).to_smiles
          end
          smi_file = Tempfile.open(['pc_ambit', '.csv'])
          pc_descriptors = nil

          # Create Ambit dataset
          smi_file.puts( "SMILES\n" )
          smi_file.puts( smiles.join("\n") )
          smi_file.close
          ambit_ds_uri = OpenTox::RestClientWrapper.post(ambit_ds_service_uri, {:file => File.new(smi_file.path)}, {:content_type => "multipart/form-data", :accept => "text/uri-list"} )
        rescue Exception => e
          LOGGER.debug "#{e.class}: #{e.message}"
          LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
        ensure
          smi_file.unlink
        end
        ambit_smiles_uri = OpenTox::RestClientWrapper.get(ambit_ds_uri + "/features", {:accept=> "text/uri-list"} ).chomp

        # Calculate 3D for CPSA
        if types.include? "cpsa"
          ambit_ds_mopac_uri = OpenTox::RestClientWrapper.post(ambit_mopac_model_uri, {:dataset_uri => ambit_ds_uri}, {:accept => "text/uri-list"} ) 
          LOGGER.debug "MOPAC dataset: #{ambit_ds_mopac_uri }"
        end

        # Get Ambit results
        ambit_result_uri = ambit_ds_uri + "?"
        ambit_result_uri << ("feature_uris[]=" + URI.encode_www_form_component(ambit_smiles_uri) + "&")
        descs_uris.each_with_index do |uri, i|
          algorithm = OpenTox::Algorithm::Generic.new(uri)
          result_uri = algorithm.run({:dataset_uri => ambit_ds_uri})
          ambit_result_uri << result_uri.split("?")[1] + "&"
          LOGGER.debug "Ambit (#{descs_uris.size}): #{i+1}"
        end
        #ambit_result_uri << ("feature_uris[]=" + URI.encode_www_form_component(ambit_smiles_uri))
        #ambit_result_uri = ("&feature_uris[]=" + URI.encode_www_form_component(ambit_smiles_uri)) + ambit_result_uri
        LOGGER.debug "Ambit result: #{ambit_result_uri}"
        ambit_result_uri

      rescue Exception => e
        LOGGER.debug "#{e.class}: #{e.message}"
        LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
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
    
    # neighbors

    module Neighbors

      # Calculate the propositionalization matrix (aka instantiation matrix) via fingerprints.
      # Same for the vector describing the query compound.
      # @param[Hash] Required keys: :neighbors, :compound, :features, :nr_hits, :fingerprints, :p_values

      def self.get_props_fingerprints (params)
        matrix = []
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
          LOGGER.debug "#{e.class}: #{e.message}"
          LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
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




    # Similarity calculations
    module Similarity

      # Tanimoto similarity
      # @param [Hash] fingerprints_a Features and values of first compound
      # @param [Hash] fingerprints_b Features and values of second compound
      # @return [Float] (Weighted) tanimoto similarity
      def self.tanimoto(fingerprints_a,fingerprints_b,weights=nil,params=nil)
        common_features = fingerprints_a.keys & fingerprints_b.keys
        all_features = (fingerprints_a.keys + fingerprints_b.keys).uniq
        if common_features.size > 0
          common_p_sum = 0.0
          common_features.each{|f| common_p_sum += [fingerprints_a[f],fingerprints_b[f]].compact.min}
          all_p_sum = 0.0
          all_features.each{|f| all_p_sum += [fingerprints_a[f],fingerprints_b[f]].compact.max}
          common_p_sum/all_p_sum
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
              dist_sum += ( (properties_a[p] - properties_b[p]) * weights[p] )**2
            else
              dist_sum += (properties_a[p] - properties_b[p])**2
            end
          end
          1/(1+Math.sqrt(dist_sum))
        else
          0.0
        end
      end

      # Cosine similarity
      # @param [Hash] properties_a Properties of first compound
      # @param [Hash] properties_b Properties of second compound
      # @return [Float] (Weighted) euclidean similarity
      def self.cosine(properties_a,properties_b,weights=nil)
        
        common_properties = properties_a.keys & properties_b.keys

        if common_properties.size > 1

          a = []; b = []
          common_properties.each do |p|
            a << properties_a[p]
            b << properties_b[p]
          end

          self.cosine_num a.to_gv, b.to_gv
        else
          0.0
        end
      end

      def self.cosine_num(a, b)
        a.dot(b) / (a.norm * b.norm)
      end

    end


  end

end

