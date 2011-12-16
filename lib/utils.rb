require 'csv'

module OpenTox
  module Algorithm


    # Calculate the propositionalization matrix (aka instantiation matrix) via fingerprints.
    # Same for the vector describing the query compound.
    # @param[Hash] Required keys: :neighbors, :compound, :features, nr_hits, :fingerprints, :p_values

    def self.get_props_fingerprints (params)
      matrix = []
      begin 

        # neighbors
        params[:neighbors].each do |n|
          n = n[:compound]
          row = []
          row_good = true
          params[:features].each do |f|
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


    # Calculate the propositionalization matrix (aka instantiation matrix) via physico-chemical descriptors.
    # Same for the vector describing the query compound.
    # The third argument takes a string from {"geometrical", "topological", "electronic", "constitutional", "hybrid" } as in ambit_descriptors.yaml
    # @param[Hash] Required keys: :dataset_uri, :pc_type
    # @return[Array, Array] Props, ids of surviving training compounds

    def self.pc_descriptors(params)

      begin
        ds = OpenTox::Dataset.find(params[:dataset_uri])
        compounds = ds.compounds.collect
        ambit_result_uri = get_pc_descriptors( { :compounds => compounds, :pc_type => params[:pc_type] } )
        csv_data = CSV.parse( OpenTox::RestClientWrapper.get(ambit_result_uri2, {:accept => "text/csv"}) )

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


  end

end


