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
      RestClientWrapper.get(uri,{:accept => "text/uri-list",:subjectid => subjectid}).to_s.each_line.collect{|u| Dataset.new(u, subjectid)}
    end

    # Load YAML representation into the dataset
    # @param [String] yaml YAML representation of the dataset
    # @return [OpenTox::Dataset] Dataset object with YAML data
    def load_yaml(yaml)
      copy YAML.load(yaml)
    end

    def load_rdfxml(rdfxml)
      raise "rdfxml data is empty" if rdfxml.to_s.size==0
      file = Tempfile.new("ot-rdfxml")
      file.puts rdfxml
      file.close
      load_rdfxml_file file
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
      if (CONFIG[:yaml_hosts].include?(URI.parse(@uri).host))
        copy YAML.load(RestClientWrapper.get(@uri, {:accept => "application/x-yaml", :subjectid => subjectid}))
      else
        parser = Parser::Owl::Dataset.new(@uri, subjectid)
        copy parser.load_uri(subjectid)
      end
    end

    # Load and return only compound URIs from the dataset service
    # @return [Array]  Compound URIs in the dataset
    def load_compounds(subjectid=nil)
      RestClientWrapper.get(File.join(uri,"compounds"),{:accept=> "text/uri-list", :subjectid => subjectid}).to_s.each_line do |compound_uri|
        @compounds << compound_uri.chomp
      end
      @compounds.uniq!
    end

    # Load and return only features from the dataset service
    # @return [Hash]  Features of the dataset
    def load_features(subjectid=nil)
      parser = Parser::Owl::Dataset.new(@uri, subjectid)
      @features = parser.load_features(subjectid)
      @features
    end

    def feature_classes(feature, subjectid=nil)
      if Feature.find(feature, subjectid).feature_type == "classification"
        classes = []
        @data_entries.each do |c,e|
          e[feature].each { |v| classes << v.to_s }
        end
        classes.uniq.sort
      else
        nil
      end
    end

=begin
    # Detect feature type(s) in the dataset
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
=end

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
    #   dataset.add "http://webservices.in-silico.ch/compound/InChI=1S/C6Cl6/c7-1-2(8)4(10)6(12)5(11)3(1)9", "http://webservices.in-silico.ch/dataset/1/feature/hamster_carcinogenicity", true
    # @param [String] compound Compound URI
    # @param [String] feature Compound URI
    # @param [Boolean,Float] value Feature value
    def add (compound,feature,value)
      @compounds << compound unless @compounds.include? compound
      @features[feature] = {}  unless @features[feature]
      @data_entries[compound] = {} unless @data_entries[compound]
      @data_entries[compound][feature] = [] unless @data_entries[compound][feature]
      @data_entries[compound][feature] << value if value!=nil
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

    # Add/modify metadata for a feature
    # @param [String] feature Feature URI
    # @param [Hash] metadata Hash with feature metadata
    def add_feature_metadata(feature,metadata)
      metadata.each { |k,v| @features[feature][k] = v }
    end
    
    # Add a new compound
    # @param [String] compound Compound URI
    def add_compound (compound)
      @compounds << compound unless @compounds.include? compound
    end
    
    # Creates a new dataset, by splitting the current dataset, i.e. using only a subset of compounds and features
    # @param [Array] compounds List of compound URIs
    # @param [Array] features List of feature URIs
    # @param [Hash] metadata Hash containing the metadata for the new dataset
    # @param [String] subjectid
    # @return [OpenTox::Dataset] newly created dataset, already saved
    def split( compounds, features, metadata, subjectid=nil)
      LOGGER.debug "split dataset using "+compounds.size.to_s+"/"+@compounds.size.to_s+" compounds"
      raise "no new compounds selected" unless compounds and compounds.size>0
      dataset = OpenTox::Dataset.create(CONFIG[:services]["opentox-dataset"],subjectid)
      if features.size==0
        compounds.each{ |c| dataset.add_compound(c) }
      else
        compounds.each do |c|
          features.each do |f|
            unless @data_entries[c][f]
              dataset.add(c,f,nil)
            else
              @data_entries[c][f].each do |v|
                dataset.add(c,f,v)
              end
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
      # TODO: rewrite feature URI's ??
      @compounds.uniq!
      if @uri
        if (CONFIG[:yaml_hosts].include?(URI.parse(@uri).host))
          RestClientWrapper.post(@uri,self.to_yaml,{:content_type =>  "application/x-yaml", :subjectid => subjectid})
        else
          File.open("ot-post-file.rdf","w+") { |f| f.write(self.to_rdfxml); @path = f.path }
          task_uri = RestClient.post(@uri, {:file => File.new(@path)},{:accept => "text/uri-list" , :subjectid => subjectid}).to_s.chomp
          #task_uri = `curl -X POST -H "Accept:text/uri-list" -F "file=@#{@path};type=application/rdf+xml" http://apps.ideaconsult.net:8080/ambit2/dataset`
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
      @data_entries[compound.uri].collect{|f,v| v.first if f.match(/prediction/)}.compact.first
    end

    def confidence(compound)
      feature_uri = @data_entries[compound.uri].collect{|f,v| f if f.match(/prediction/)}.compact.first
      @features[feature_uri][OT.confidence]
    end

    def descriptors(compound)
      @data_entries[compound.uri].collect{|f,v| @features[f] if f.match(/descriptor/)}.compact if @data_entries[compound.uri]
    end

    def measured_activities(compound)
      source = @metadata[OT.hasSource]
      @data_entries[compound.uri].collect{|f,v| v if f.match(/#{source}/)}.compact.flatten
    end

    def neighbors(compound)
      @data_entries[compound.uri].collect{|f,v| @features[f] if f.match(/neighbor/)}.compact
    end

#    def errors(compound)
#      features = @data_entries[compound.uri].keys
#      features.collect{|f| @features[f][OT.error]}.join(" ") if features
#    end

  end
end
