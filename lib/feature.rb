module OpenTox
  class Feature
    include OpenTox

    attr_accessor :subjectid

    # Find a feature
    # @param [String] uri Feature URI
    # @return [OpenTox::Task] Feature object
    def self.find(uri, subjectid=nil)
      return nil unless uri   
      feature = Feature.new uri
      if (CONFIG[:yaml_hosts].include?(URI.parse(uri).host))
        feature.add_metadata YAML.load(RestClientWrapper.get(uri,{:accept => "application/x-yaml", :subjectid => subjectid}))
      else
        feature.add_metadata  Parser::Owl::Dataset.new(uri).load_metadata
      end
      feature.subjectid = subjectid
      feature
    end

    # provides feature type, possible types are "regression" or "classification"
    # @return [String] feature type, unknown if OT.isA property is unknown/ not set
    def feature_type
      raise OpenTox::BadRequestError.new("rdf type of feature '"+uri.to_s+"' not set") unless metadata[RDF.type]
      if metadata[RDF.type].to_a.flatten.include?(OT.NominalFeature)
        "classification"
      elsif metadata[RDF.type].to_a.flatten.include?(OT.NumericFeature)
        "regression"
      elsif metadata[OWL.sameAs]
        metadata[OWL.sameAs].each do |f|
          begin
            type = Feature.find(f, subjectid).feature_type
            return type unless type=="unknown"
          rescue => ex
            LOGGER.warn "could not load same-as-feature '"+f.to_s+"' for feature '"+uri.to_s+"' : "+ex.message.to_s
          end
        end
        "unknown"
      else
        "unknown"
      end
    end    
  end

  # Get OWL-DL representation in RDF/XML format
  # @return [application/rdf+xml] RDF/XML representation
  def to_rdfxml
    s = Serializer::Owl.new
    s.add_feature(@uri,@metadata)
    @metadata.values.grep(/model\/\d+$/).each{ |m| s.add_uri(m,OT.Model)}
    @metadata.values.grep(/feature/).each{ |f| s.add_uri(f,OT.Feature)}
    #s.add_parameters(@uri,@parameters) if @parameters
    s.to_rdfxml
  end
end
