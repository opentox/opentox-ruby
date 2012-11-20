module OpenTox
  
  # Ruby wrapper for OpenTox Dataset Webservices (http://opentox.org/dev/apis/api-1.2/dataset).
  class Dataset 

    include OpenTox

    attr_reader :features, :compounds, :data_entries, :metadata

    # Create dataset with optional URI. Does not load data into the dataset - you will need to execute one of the load_* methods to pull data from a service or to insert it from other representations.
    # @example Create an empty dataset
    #   dataset = OpenTox::Dataset.new
    # @example Create an empty dataset with URI
    #   dataset = OpenTox::Dataset.new("http:://webservices.in-silico/ch/dataset/1")
    # @param [optional, String] uri Dataset URI
    # @return [OpenTox::Dataset] Dataset object
    def initialize(uri=nil,subjectid=nil)
      super uri
      @features = {}
      @compounds = []
      @data_entries = {}
    end

    # Create an empty dataset and save it at the dataset service (assigns URI to dataset)
    # @example Create new dataset and save it to obtain a URI 
    #   dataset = OpenTox::Dataset.create
    # @param [optional, String] uri Dataset URI
    # @return [OpenTox::Dataset] Dataset object
    def self.create(uri=CONFIG[:services]["opentox-dataset"], subjectid=nil)
      dataset = Dataset.new(nil,subjectid)
      dataset.save(subjectid)
      dataset
    end

    # Create dataset from CSV file (format specification: http://toxcreate.org/help)
    # - loads data_entries, compounds, features
    # - sets metadata (warnings) for parser errors
    # - you will have to set remaining metadata manually
    # @param [String] file CSV file path
    # @return [OpenTox::Dataset] Dataset object with CSV data
    def self.create_from_csv_file(file, subjectid=nil) 
      dataset = Dataset.create(CONFIG[:services]["opentox-dataset"], subjectid)
      parser = Parser::Spreadsheets.new
      parser.dataset = dataset
      parser.load_csv(File.open(file).read)
      dataset.save(subjectid)
      dataset
    end

    def self.from_json(json, subjectid=nil) 
      dataset = Dataset.new(nil,subjectid)
      dataset.copy_hash Yajl::Parser.parse(json)
      dataset
    end
    
    # Find a dataset and load all data. This can be time consuming, use Dataset.new together with one of the load_* methods for a fine grained control over data loading.
    # @param [String] uri Dataset URI
    # @return [OpenTox::Dataset] Dataset object with all data
    def self.find(uri, subjectid=nil)
      return nil unless uri
      dataset = Dataset.new(uri, subjectid)
      dataset.load_all(subjectid)
      dataset
    end
    
    # replaces find as exist check, takes not as long, does NOT raise an un-authorized exception
    # @param [String] uri Dataset URI
    # @return [Boolean] true if dataset exists and user has get rights, false else 
    def self.exist?(uri, subjectid=nil)
      return false unless uri
      dataset = Dataset.new(uri, subjectid)
      begin
        dataset.load_metadata( subjectid ).size > 0
      rescue
        false
      end
    end

    # Get all datasets from a service
    # @param [optional,String] uri URI of the dataset service, defaults to service specified in configuration
    # @return [Array] Array of dataset object without data (use one of the load_* methods to pull data from the server)
    def self.all(uri=CONFIG[:services]["opentox-dataset"], subjectid=nil)
      RestClientWrapper.get(uri,{:accept => "text/uri-list",:subjectid => subjectid}).to_s.each_line.collect{|u| Dataset.new(u.chomp, subjectid)}
    end

    # Load YAML representation into the dataset
    # @param [String] yaml YAML representation of the dataset
    # @return [OpenTox::Dataset] Dataset object with YAML data
    def load_yaml(yaml)
      copy YAML.load(yaml)
    end

    def load_json(json)
      copy_hash Yajl::Parser.parse(json)
    end

    def load_rdfxml(rdfxml, subjectid=nil)
      raise "rdfxml data is empty" if rdfxml.to_s.size==0
      file = Tempfile.new("ot-rdfxml")
      file.puts rdfxml
      file.close
      load_rdfxml_file file, subjectid
      file.delete
    end

    # Load RDF/XML representation from a file
    # @param [String] file File with RDF/XML representation of the dataset
    # @return [OpenTox::Dataset] Dataset object with RDF/XML data
    def load_rdfxml_file(file, subjectid=nil)
      parser = Parser::Owl::Dataset.new @uri, subjectid
      parser.uri = file.path
      copy parser.load_uri(subjectid)
    end

    def load_sdf(sdf,subjectid=nil)
      save(subjectid) unless @uri # get a uri for creating features
      parser = Parser::Sdf.new
      parser.dataset = self
      parser.load_sdf(sdf)
    end

    # Load CSV string (format specification: http://toxcreate.org/help)
    # - loads data_entries, compounds, features
    # - sets metadata (warnings) for parser errors
    # - you will have to set remaining metadata manually
    # @param [String] csv CSV representation of the dataset
    # @return [OpenTox::Dataset] Dataset object with CSV data
    def load_csv(csv, subjectid=nil) 
      save(subjectid) unless @uri # get a uri for creating features
      parser = Parser::Spreadsheets.new
      parser.dataset = self
      parser.load_csv(csv)
    end

    # Load Spreadsheet book (created with roo gem http://roo.rubyforge.org/, excel format specification: http://toxcreate.org/help)
    # - loads data_entries, compounds, features
    # - sets metadata (warnings) for parser errors
    # - you will have to set remaining metadata manually
    # @param [Excel] book Excel workbook object (created with roo gem)
    # @return [OpenTox::Dataset] Dataset object with Excel data
    def load_spreadsheet(book, subjectid=nil)
      save(subjectid) unless @uri # get a uri for creating features
      parser = Parser::Spreadsheets.new
      parser.dataset = self
      parser.load_spreadsheet(book)
    end
    
    # Load and return only metadata of a Dataset object
    # @return [Hash] Metadata of the dataset
    def load_metadata(subjectid=nil)
      add_metadata Parser::Owl::Dataset.new(@uri, subjectid).load_metadata(subjectid)
      self.uri = @uri if @uri # keep uri
      @metadata
    end

    # Load all data (metadata, data_entries, compounds and features) from URI
    def load_all(subjectid=nil)
      if (CONFIG[:json_hosts].include?(URI.parse(@uri).host))
        copy_hash Yajl::Parser.parse(RestClientWrapper.get(@uri, {:accept => "application/json", :subjectid => subjectid}))
      else
        parser = Parser::Owl::Dataset.new(@uri, subjectid)
        copy parser.load_uri(subjectid)
      end
    end

    # Load and return only compound URIs from the dataset service
    # @return [Array]  Compound URIs in the dataset
    def load_compounds(subjectid=nil)
      # fix for datasets like http://apps.ideaconsult.net:8080/ambit2/dataset/272?max=50
      u = URI::parse(uri)
      u.path = File.join(u.path,"compounds")
      u = u.to_s
      RestClientWrapper.get(u,{:accept=> "text/uri-list", :subjectid => subjectid}).to_s.each_line do |compound_uri|
        @compounds << compound_uri.chomp
      end
      @compounds
    end

    # Load and return only features from the dataset service
    # @return [Hash]  Features of the dataset
    def load_features(subjectid=nil)
      if (CONFIG[:json_hosts].include?(URI.parse(@uri).host))
        @features = Yajl::Parser.parse(RestClientWrapper.get(File.join(@uri,"features"), {:accept => "application/json", :subjectid => subjectid}))
      else
        parser = Parser::Owl::Dataset.new(@uri, subjectid)
        @features = parser.load_features(subjectid)
      end
      @features
    end

    # returns the accept_values of a feature, i.e. the classification domain / all possible feature values 
    # @param [String] feature the URI of the feature
    # @return [Array] return array with strings, nil if value is not set (e.g. when feature is numeric)
    def accept_values(feature)
      accept_values = features[feature][OT.acceptValue]
      if accept_values
        accept_values = accept_values.collect{|av| av.to_s}
        accept_values.sort!
      end 
      accept_values
    end

    # Detect feature type (reduced to one across all features)
    # Classification takes precedence over regression
    # DEPRECATED -- 
    #   HAS NO SENSE FOR DATASETS WITH MORE THAN 1 FEATURE
    #   FEATURES CAN HAVE MULTIPLE TYPES
    # Replacement: see feature_types()
    # @return [String] `classification", "regression", "mixed" or unknown`
    def feature_type(subjectid=nil)
      load_features(subjectid)
      feature_types = @features.collect{|f,metadata| metadata[RDF.type]}.flatten.uniq
      if feature_types.include?(OT.NominalFeature)
        "classification"
      elsif feature_types.include?(OT.NumericFeature)
        "regression"
      else
        "unknown"
      end
    end


    # Detect feature types. A feature can have multiple types.
    # Returns types hashed by feature URI, with missing features omitted.
    # Example (YAML):
    #   http://toxcreate3.in-silico.ch:8082/dataset/152/feature/nHal: 
    #   - http://www.opentox.org/api/1.1#NumericFeature
    #   - http://www.opentox.org/api/1.1#NominalFeature
    #   ...
    #
    # @return [Hash] Keys: feature URIs, Values: Array of types
    def feature_types(subjectid=nil)
      load_features(subjectid)
        @features.inject({}){ |h,(f,metadata)| 
          h[f]=metadata[RDF.type] unless metadata[RDF.type][0].include? "MissingFeature" 
          h
        }
    end
=begin
=end

    def to_json
      Yajl::Encoder.encode({:uri => @uri, :metadata => @metadata, :data_entries => @data_entries, :compounds => @compounds, :features => @features})
    end

    # Get Spreadsheet representation
    # @return [Spreadsheet::Workbook] Workbook which can be written with the spreadsheet gem (data_entries only, metadata will will be discarded))
    def to_spreadsheet
      Serializer::Spreadsheets.new(self).to_spreadsheet
    end

    # Get Excel representation (alias for to_spreadsheet)
    # @return [Spreadsheet::Workbook] Workbook which can be written with the spreadsheet gem (data_entries only, metadata will will be discarded))
    def to_xls
      to_spreadsheet
    end

    # Get CSV string representation (data_entries only, metadata will be discarded)
    # @return [String] CSV representation
    def to_csv
      Serializer::Spreadsheets.new(self).to_csv
    end

    # Get OWL-DL in ntriples format
    # @return [String] N-Triples representation
    def to_ntriples
      s = Serializer::Owl.new
      s.add_dataset(self)
      s.to_ntriples
    end

    # Get OWL-DL in RDF/XML format
    # @return [String] RDF/XML representation
    def to_rdfxml
      s = Serializer::Owl.new
      s.add_dataset(self)
      s.to_rdfxml
    end

    # Get SDF representation of compounds
    # @return [String] SDF representation
    def to_sdf
      sum=""
      @compounds.each{ |c|
        sum << OpenTox::Compound.new(c).to_inchi
        sum << OpenTox::Compound.new(c).to_sdf.sub(/\n\$\$\$\$/,'')
        @data_entries[c].each{ |f,v|
          sum << ">  <\"#{f}\">\n"
          sum << v.join(", ")
          sum << "\n\n"
        }
        sum << "$$$$\n"
      }
      sum
    end

    def to_urilist
      @compounds.inject { |sum, c|
        sum << OpenTox::Compound.new(c).uri
        sum + "\n"
      }
    end

    # Get name (DC.title) of a feature
    # @param [String] feature Feature URI
    # @return [String] Feture title
    def feature_name(feature)
      @features[feature][DC.title]
    end

    def title
      @metadata[DC.title]
    end

    # Insert a statement (compound_uri,feature_uri,value)
    # @example Insert a statement (compound_uri,feature_uri,value)
    #   dataset.add "http://webservices.in-silico.ch/compound/InChI=1S/C6Cl6/c7-1-2(8)4(10)6(12)5(11)3(1)9", "http://webservices.in-silico.ch/dataset/1/feature/hamster_carcinogenicity", 1
    # @param [String] compound Compound URI
    # @param [String] feature Compound URI
    # @param [Boolean,Float] value Feature value
    def add (compound,feature,value)
      LOGGER.warn("dataset.add is deprecated and should not be used any longer")
      self.add_compound(compound)
      self.add_data_entry(compound,feature,value)
    end

    # Insert a data entry
    # @param [String] compound Compound URI
    # @param [String] feature Compound URI
    # @param [Boolean,Float] value Feature value
    def add_data_entry (compound,feature,value)
      if @compounds.include? compound
        @features[feature] = {}  unless @features[feature]
        @data_entries[compound] = {} unless @data_entries[compound]
        @data_entries[compound][feature] = [] unless @data_entries[compound][feature]
        @data_entries[compound][feature] << value
      end
    end


    # Add/modify metadata, existing entries will be overwritten
    # @example
    #   dataset.add_metadata({DC.title => "any_title", DC.creator => "my_email"})
    # @param [Hash] metadata Hash mapping predicate_uris to values
    def add_metadata(metadata)
      metadata.each { |k,v| @metadata[k] = v }
    end

    # Add a feature
    # @param [String] feature Feature URI
    # @param [Hash] metadata Hash with feature metadata
    def add_feature(feature,metadata={})
      @features[feature] = metadata
    end

    # Complete feature values by adding zeroes
    # @param [Hash] key: compound, value: duplicate sizes
    def complete_data_entries(compound_sizes)
      all_features = @features.keys
      @data_entries.each { |c, e|
        (Set.new(all_features.collect)).subtract(Set.new e.keys).to_a.each { |f|
          compound_sizes[c].times { 
            self.add(c,f,0) 
          }
        }
      }
    end

    # Add/modify metadata for a feature
    # @param [String] feature Feature URI
    # @param [Hash] metadata Hash with feature metadata
    def add_feature_metadata(feature,metadata)
      metadata.each { |k,v| @features[feature][k] = v }
    end
    
    # Add a new compound
    # @param [String] compound Compound URI
    def add_compound (compound)
      @compounds << compound
      #@compounds << compound unless @compounds.include? compound
    end
    
    protected
    # returns array with (row-)indices of a compound (in the compounds arrays) 
    def compound_indices( compound )
      unless defined?(@cmp_indices) and @cmp_indices.has_key?(compound)
        @cmp_indices = {}
        @compounds.size.times do |i|
          c = @compounds[i]
          if @cmp_indices[c]==nil
            @cmp_indices[c] = [i]
          else
            @cmp_indices[c] = @cmp_indices[c]+[i]
          end   
        end
      end
      @cmp_indices[compound]
    end
    
    # returns index in data_entries-array for the compound-index 
    def entry_index( compound_index )
      unless defined?(@entry_indices) and @entry_indices.has_key?(compound_index)
        @entry_indices = {}
        @compounds.size.times do |i|
          @entry_indices[i] = compound_indices(@compounds[i]).index(i)
        end
      end
      @entry_indices[compound_index]
    end
    public 

    # maps a compound-index from another dataset to a compound-index from this dataset
    # mapping works as follows:
    # (compound c is the compound identified by the compound-index of the other dataset)
    # * c occurs only once in this dataset? map compound-index of other dataset to index in this dataset
    # * c occurs >1 in this dataset?
    # ** number of occurences is equal in both datasets? assume order is preserved(!) and map accordingly
    # ** number of occurences is not equal in both datasets? cannot map, raise error
    # @param [OpenTox::Dataset] dataset that should be mapped to this dataset (fully loaded)
    # @param [Fixnum] compound_index, corresponding to dataset
    def compound_index( dataset, compound_index )
      unless defined?(@index_map) and @index_map[dataset.uri]
        map = {}
        dataset.compounds.uniq.each do |compound|
          self_indices = compound_indices(compound)
          next unless self_indices
          dataset_indices = dataset.compound_indices(compound)
          if self_indices.size==1  
            dataset_indices.size.times do |i|
              map[dataset_indices[i]] = self_indices[0]
            end
          elsif self_indices.size==dataset_indices.size
            # we do assume that the order is preseverd!
            dataset_indices.size.times do |i|
              map[dataset_indices[i]] = self_indices[i]
            end
          else
            raise "cannot map compound #{compound} from dataset #{dataset.uri} to dataset #{uri}, "+
              "compound occurs #{dataset_indices.size} times and #{self_indices.size} times"
          end
        end  
        @index_map = {} unless defined?(@index_map)
        @index_map[dataset.uri] = map
      end
      @index_map[dataset.uri][compound_index]
    end
    
    # returns data entry value for a compound-index and given feature
    # @param [Fixnum] compound_index
    # @param [OpenTox::Feature] feature
    def data_entry_value( compound_index, feature )
      raise "please give compound index instead of '#{compound_index}'" unless compound_index.is_a?(Fixnum)
      c = @compounds[compound_index]
      if @data_entries[c]==nil
        nil
      elsif @data_entries[c][feature]==nil
        nil
      else
        @data_entries[c][feature][entry_index(compound_index)]
      end  
    end
    
    # Creates a new dataset, by splitting the current dataset, i.e. using only a subset of compounds and features
    # @param [Array] compound_indices List of compound indices
    # @param [Array] features List of feature URIs
    # @param [Hash] metadata Hash containing the metadata for the new dataset
    # @param [String] subjectid
    # @return [OpenTox::Dataset] newly created dataset, already saved
    def split( compound_indices, features, metadata, subjectid=nil)
      raise "Dataset.split : pls give compounds as indices" if compound_indices.size==nil or !compound_indices[0].is_a?(Fixnum) 
      LOGGER.debug "split dataset using "+compound_indices.size.to_s+"/"+@compounds.size.to_s+" compounds"
      dataset = OpenTox::Dataset.create(CONFIG[:services]["opentox-dataset"],subjectid)
      features.each{|f| dataset.add_feature(f,@features[f])}
      compound_indices.each do |c_idx|
        c = @compounds[c_idx]
        dataset.add_compound(c)
        if @data_entries[c]
          features.each do |f|
            if @data_entries[c][f] 
              dataset.add_data_entry c,f,@data_entries[c][f][entry_index(c_idx)]
            else
              dataset.add_data_entry c,f,nil
            end
          end
        end
      end
      dataset.add_metadata(metadata)
      dataset.save(subjectid)
      dataset
    end
    
    # Save dataset at the dataset service 
    # - creates a new dataset if uri is not set
    # - overwrites dataset if uri exists
    # @return [String] Dataset URI
    def save(subjectid=nil)
      if @uri
        if (CONFIG[:json_hosts].include?(URI.parse(@uri).host))
          #LOGGER.debug self.to_json
          RestClientWrapper.post(@uri,self.to_json,{:content_type =>  "application/json", :subjectid => subjectid})
        else
          File.open("ot-post-file.rdf","w+") { |f| f.write(self.to_rdfxml); @path = f.path }
          task_uri = RestClient.post(@uri, {:file => File.new(@path)},{:accept => "text/uri-list" , :subjectid => subjectid}).to_s.chomp
          Task.find(task_uri).wait_for_completion
          self.uri = RestClientWrapper.get(task_uri,{:accept => 'text/uri-list', :subjectid => subjectid})
        end
      else
        # create dataset if uri is empty
        self.uri = RestClientWrapper.post(CONFIG[:services]["opentox-dataset"],{:subjectid => subjectid}).to_s.chomp
      end
      @uri
    end

    # Delete dataset at the dataset service
    def delete(subjectid=nil)
      RestClientWrapper.delete(@uri, :subjectid => subjectid)
    end

    # Copy a hash (eg. from JSON) into a dataset (rewrites URI)
    def copy_hash(hash)
      @metadata = hash["metadata"]
      @data_entries = hash["data_entries"]
      @compounds = hash["compounds"]
      @features = hash["features"]
      if @uri
        self.uri = @uri 
      else
        @uri = hash["metadata"][XSD.anyURI]
      end
    end

    def value_map(prediction_feature_uri)
      training_classes = accept_values(prediction_feature_uri).sort
      value_map=Hash.new 
      training_classes.each_with_index { |c,i| value_map[i+1] = c }
      value_map
    end


    private
    # Copy a dataset (rewrites URI)
    def copy(dataset)
      @metadata = dataset.metadata
      @data_entries = dataset.data_entries
      @compounds = dataset.compounds
      @features = dataset.features
      if @uri
        self.uri = @uri 
      else
        @uri = dataset.metadata[XSD.anyURI]
      end
    end
  end

  # Class with special methods for lazar prediction datasets
  class LazarPrediction < Dataset

    # Find a prediction dataset and load all data. 
    # @param [String] uri Prediction dataset URI
    # @return [OpenTox::Dataset] Prediction dataset object with all data
    def self.find(uri, subjectid=nil)
      prediction = LazarPrediction.new(uri, subjectid)
      prediction.load_all(subjectid)
      prediction
    end

    def value(compound)
      v = nil
      v = @data_entries[compound.uri].collect{|f,v| v.first if f.match(/value/)}.compact.first if @data_entries[compound.uri]
      v = nil if v.is_a? Array and v.empty?
      v
    end

    def confidence(compound)
      @data_entries[compound.uri].collect{|f,v| v.first if f.match(/confidence/)}.compact.first if @data_entries[compound.uri]
    end

    def descriptors(compound)
      @data_entries[compound.uri].collect{|f,v| @features[f] if f.match(/descriptor/)}.compact if @data_entries[compound.uri]
    end

    def measured_activities(compound)
      @data_entries[compound.uri].collect{|f,v| v if f.match(/#{@metadata[OT.hasSource]}/)}.compact.flatten if @data_entries[compound.uri]
    end

    def neighbors(compound)
      @data_entries[compound.uri].collect{|f,v| @features[f] if f.match(/neighbor/)}.compact if @data_entries[compound.uri]
    end


#    def errors(compound)
#      features = @data_entries[compound.uri].keys
#      features.collect{|f| @features[f][OT.error]}.join(" ") if features
#    end

  end
end
