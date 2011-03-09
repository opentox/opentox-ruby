module OpenTox

  attr_reader :uri
  attr_accessor :metadata

  # Initialize OpenTox object with optional uri
  # @param [optional, String] URI
  def initialize(uri=nil)
    @metadata = {}
    self.uri = uri if uri
  end

  # Set URI
  # @param [String] URI
  def uri=(uri)
    @uri = uri
    @metadata[XSD.anyURI] = uri
  end

  # Get all objects from a service
  # @return [Array] List of available URIs
  def self.all(uri, subjectid=nil)
    RestClientWrapper.get(uri,:accept => "text/uri-list", :subjectid => subjectid).to_s.split(/\n/)
  end

  # Load (and return) metadata from object URI
  # @return [Hash] Metadata
  def load_metadata(subjectid=nil)
    @metadata = Parser::Owl::Generic.new(@uri).load_metadata(subjectid)
    @metadata
  end

  def add_metadata(metadata)
    metadata.each { |k,v| @metadata[k] = v }
  end

  # Get OWL-DL representation in RDF/XML format
  # @return [application/rdf+xml] RDF/XML representation
  def to_rdfxml
    s = Serializer::Owl.new
    s.add_metadata(@uri,@metadata)
    #s.add_parameters(@uri,@parameters) if @parameters
    s.to_rdfxml
  end

  # deletes the resource, deletion should have worked when no RestCallError raised
  def delete(subjectid=nil)
    RestClientWrapper.delete(uri,:subjectid => subjectid)
  end

end

