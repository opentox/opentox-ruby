require 'csv'


module OpenTox

  module Algorithm

    include OpenTox

    # Calculate physico-chemical descriptors.
    # @param[Hash] Required keys: :dataset_uri, :pc_type
    # @return[String] dataset uri

    def self.pc_descriptors(params)

      begin
        ds = OpenTox::Dataset.find(params[:dataset_uri])
        compounds = ds.compounds.collect
        ambit_result_uri, smiles_to_inchi = get_pc_descriptors( { :compounds => compounds, :pc_type => params[:pc_type] } )
        #ambit_result_uri = ["http://apps.ideaconsult.net:8080/ambit2/dataset/987103?" ,"feature_uris[]=http%3A%2F%2Fapps.ideaconsult.net%3A8080%2Fambit2%2Ffeature%2F4276789&", "feature_uris[]=http%3A%2F%2Fapps.ideaconsult.net%3A8080%2Fambit2%2Fmodel%2F16%2Fpredicted"] # for testing
        LOGGER.debug "Ambit result uri for #{params.inspect}: '#{ambit_result_uri.to_yaml}'"
        load_ds_csv(ambit_result_uri, smiles_to_inchi)
      rescue Exception => e
        LOGGER.debug "#{e.class}: #{e.message}"
        LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
      end

    end
    
    # Calculates PC descriptors via Ambit -- DO NOT OVERLOAD Ambit.
    # @param[Hash] Required keys: :compounds, :pc_type
    # @return[Array] Ambit result uri, piecewise (1st: base, 2nd: SMILES, 3rd+: features
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
          smiles_array = []; smiles_to_inchi = {}
          params[:compounds].each do |n|
            cmpd = OpenTox::Compound.new(n)
            smiles_string = cmpd.to_smiles
            smiles_to_inchi[smiles_string] = URI.encode_www_form_component(cmpd.to_inchi)
            smiles_array << smiles_string
          end
          smi_file = Tempfile.open(['pc_ambit', '.csv'])
          pc_descriptors = nil

          # Create Ambit dataset
          smi_file.puts( "SMILES\n" )
          smi_file.puts( smiles_array.join("\n") )
          smi_file.flush
          ambit_ds_uri = OpenTox::RestClientWrapper.post(ambit_ds_service_uri, {:file => File.new(smi_file.path)}, {:content_type => "multipart/form-data", :accept => "text/uri-list"} )
        rescue Exception => e
          LOGGER.debug "#{e.class}: #{e.message}"
          LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
        ensure
          smi_file.close! if smi_file
        end
        ambit_smiles_uri = OpenTox::RestClientWrapper.get(ambit_ds_uri + "/features", {:accept=> "text/uri-list"} ).chomp

        # Calculate 3D for CPSA
        if types.include? "cpsa"
          ambit_ds_mopac_uri = OpenTox::RestClientWrapper.post(ambit_mopac_model_uri, {:dataset_uri => ambit_ds_uri}, {:accept => "text/uri-list"} ) 
          LOGGER.debug "MOPAC dataset: #{ambit_ds_mopac_uri }"
        end

        # Get Ambit results
        ambit_result_uri = [] # 1st pos: base uri, then features
        ambit_result_uri << ambit_ds_uri + "?"
        ambit_result_uri << ("feature_uris[]=" + URI.encode_www_form_component(ambit_smiles_uri) + "&")
        descs_uris.each_with_index do |uri, i|
          algorithm = Algorithm::Generic.new(uri)
          result_uri = algorithm.run({:dataset_uri => ambit_ds_uri})
          ambit_result_uri << result_uri.split("?")[1] + "&"
          LOGGER.debug "Ambit (#{descs_uris.size}): #{i+1}"
        end
        #LOGGER.debug "Ambit result: #{ambit_result_uri.join('')}"
        [ ambit_result_uri, smiles_to_inchi ]

      rescue Exception => e
        LOGGER.debug "#{e.class}: #{e.message}"
        LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
      end
    end


    # Load dataset via CSV
    # @param[Array] Ambit result uri, piecewise (1st: base, 2nd: SMILES, 3rd+: features
    # @return[String] dataset uri 
    def self.load_ds_csv(ambit_result_uri, smiles_to_inchi)
      
      master=nil
      (1...ambit_result_uri.size).collect { |idx|
        curr_uri = ambit_result_uri[0] + ambit_result_uri[idx]
        LOGGER.debug "Requesting #{curr_uri}"
        csv_data = CSV.parse( OpenTox::RestClientWrapper.get(curr_uri, {:accept => "text/csv"}) )
        if csv_data[0] && csv_data[0].size>1
          if master.nil? # This is the smiles entry
            (1...csv_data.size).each{ |idx| csv_data[idx][1] = smiles_to_inchi[csv_data[idx][1]] }
            master = csv_data
            next
          else
            index_uri = csv_data[0].index("SMILES")
            csv_data.map {|i| i.delete_at(index_uri)} if index_uri #Removes additional SMILES information
            
            nr_cols = (csv_data[0].size)-1
            LOGGER.debug "Merging #{nr_cols} new columns"
            master.each {|row| nr_cols.times { row.push(nil) }  } # Adds empty columns to all rows
            csv_data.each do |row|
              temp = master.assoc(row[0]) # Finds the appropriate line in master
              ((-1*nr_cols)..-1).collect.each { |idx|
                temp[idx] = row[nr_cols+idx+1] if temp # Updates columns if line is found
              }
            end
          end
        end
      }

      index_uri = master[0].index("Compound")
      master.map {|i| i.delete_at(index_uri)}
      master[0].each {|cell| cell.chomp!(" ")}
      master[0][0] = "Compound" #"SMILES" 
      index_smi = master[0].index("SMILES")
      master.map {|i| i.delete_at(index_smi)} if index_smi
      #master[0][0] = "SMILES" 
       
      #LOGGER.debug "-------- AM: Writing to dumpfile"
      #File.open("/tmp/test.csv", 'w') {|f| f.write( master.collect {|r| r.join(",")}.join("\n") ) }
     
      parser = OpenTox::Parser::Spreadsheets.new
      ds = OpenTox::Dataset.new
      ds.save
      parser.dataset = ds
      ds = parser.load_csv(master.collect{|r| r.join(",")}.join("\n"))
      ds.save
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
      return array.uniq.size == 1
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
      
      # Get confidence.
      # @param[Hash] Required keys: :sims, :acts
      # @return[Float] Confidence
      def self.get_confidence(params)
        conf = params[:sims].inject{|sum,x| sum + x }
        confidence = conf/params[:sims].size
        LOGGER.debug "Confidence is: '" + confidence.to_s + "'."
        return confidence
      end

    end


    # Similarity calculations
    module Similarity

      # Tanimoto similarity
      # @param [Hash, Array] fingerprints of first compound
      # @param [Hash, Array] fingerprints of second compound
      # @return [Float] (Weighted) tanimoto similarity
      def self.tanimoto(fingerprints_a,fingerprints_b,weights=nil,params=nil)

        common_p_sum = 0.0
        all_p_sum = 0.0

        # fingerprints are hashes
        if fingerprints_a.class == Hash && fingerprints_b.class == Hash
          common_features = fingerprints_a.keys & fingerprints_b.keys
          all_features = (fingerprints_a.keys + fingerprints_b.keys).uniq
          if common_features.size > 0
            common_features.each{ |f| common_p_sum += [ fingerprints_a[f], fingerprints_b[f] ].min }
            all_features.each{ |f| all_p_sum += [ fingerprints_a[f],fingerprints_b[f] ].compact.max } # compact, since one fp may be empty at that pos
          end

        # fingerprints are arrays
        elsif fingerprints_a.class == Array && fingerprints_b.class == Array
          size = [ fingerprints_a.size, fingerprints_b.size ].min
          LOGGER.warn "fingerprints don't have equal size" if fingerprints_a.size != fingerprints_b.size
          (0...size).each { |idx|
            common_p_sum += [ fingerprints_a[idx], fingerprints_b[idx] ].min
            all_p_sum += [ fingerprints_a[idx], fingerprints_b[idx] ].max
          }
        end

        (all_p_sum > 0.0) ? (common_p_sum/all_p_sum) : 0.0

      end


      # Cosine similarity
      # @param [Hash] properties_a key-value properties of first compound
      # @param [Hash] properties_b key-value properties of second compound
      # @return [Float] cosine of angle enclosed between vectors induced by keys present in both a and b
      def self.cosine(fingerprints_a,fingerprints_b,weights=nil)

        # fingerprints are hashes
        if fingerprints_a.class == Hash && fingerprints_b.class == Hash
          a = []; b = []
          common_features = fingerprints_a.keys & fingerprints_b.keys
          if common_features.size > 1
            common_features.each do |p|
              a << fingerprints_a[p]
              b << fingerprints_b[p]
            end
          end

        # fingerprints are arrays
        elsif fingerprints_a.class == Array && fingerprints_b.class == Array
          a = fingerprints_a
          b = fingerprints_b
        end

        (a.size > 0 && b.size > 0) ? self.cosine_num(a.to_gv, b.to_gv) : 0.0

      end


      # Cosine similarity
      # @param [GSL::Vector] a
      # @param [GSL::Vector] b
      # @return [Float] cosine of angle enclosed between a and b
      def self.cosine_num(a, b)
        if a.size>12 && b.size>12
          a = a[0..11]
          b = b[0..11]
        end
        a.dot(b) / (a.norm * b.norm)
      end


      # Outlier detection based on Mahalanobis distances
      # Multivariate detection on X, univariate detection on y
      # Uses an existing Rinruby instance, if possible
      # @param[Hash] Keys query_matrix, data_matrix, acts are required; r, p_outlier optional
      # @return[Array] indices identifying outliers (may occur several times, this is intended)
      def self.outliers(params)
        outlier_array = []
        data_matrix = params[:data_matrix]
        query_matrix = params[:query_matrix]
        acts = params[:acts]
        begin
          LOGGER.debug "Outliers (p=#{params[:p_outlier] || 0.9999})..."
          r = ( params[:r] || RinRuby.new(false,false) )
          r.eval "suppressPackageStartupMessages(library(\"robustbase\"))"
          r.eval "outlier_threshold = #{params[:p_outlier] || 0.999}"
          nr_cases, nr_features = data_matrix.to_a.size, data_matrix.to_a[0].size
          r.odx = data_matrix.to_a.flatten
          r.q = query_matrix.to_a.flatten
          r.y = acts.to_a.flatten
          r.eval "odx = matrix(odx, #{nr_cases}, #{nr_features}, byrow=T)"
          r.eval 'odx = rbind(q,odx)' # query is nr 0 (1) in ruby (R)
          r.eval 'mah = covMcd(odx)$mah' # run MCD alg
          r.eval "mah = pchisq(mah,#{nr_features})"
          r.eval 'outlier_array = which(mah>outlier_threshold)'  # multivariate outliers using robust mahalanobis
          outlier_array = r.outlier_array.to_a.collect{|v| v-2 }  # translate to ruby index (-1 for q, -1 due to ruby)
          r.eval 'fqu = matrix(summary(y))[2]'
          r.eval 'tqu = matrix(summary(y))[5]'
          r.eval 'outlier_array = which(y>(tqu+1.5*IQR(y)))'     # univariate outliers due to Tukey (http://goo.gl/mwzNH)
          outlier_array += r.outlier_array.to_a.collect{|v| v-1 } # translate to ruby index (-1 due to ruby)
          r.eval 'outlier_array = which(y<(fqu-1.5*IQR(y)))'
          outlier_array += r.outlier_array.to_a.collect{|v| v-1 }
        rescue Exception => e
          LOGGER.debug "#{e.class}: #{e.message}"
          #LOGGER.debug "Backtrace:\n\t#{e.backtrace.join("\n\t")}"
        end
        outlier_array
      end

    end


  end

end

