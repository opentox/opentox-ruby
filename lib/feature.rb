module OpenTox
  class Feature
    include OpenTox

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
      feature
    end
    
    # provides feature type, possible types are "regression" or "classification"
    # @return [String] feature type, unknown if OT.isA property is unknown/ not set
    def feature_type
      if metadata[RDF.type].flatten.include?(OT.NominalFeature)
        "classification"
      elsif metadata[RDF.type].flatten.include?(OT.NumericFeature)
        "regression"
      else
        #"unknown"
        metadata[RDF.type].inspect
      end
=begin
      case metadata[RDF.type]
      when /NominalFeature/
        "classification"
      when /NumericFeature/
        "regression"
      else
        "unknown"
      end
=end
    end    
    
  end
end
