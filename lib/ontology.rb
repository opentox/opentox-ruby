module OpenTox
  module Ontology
    module Echa

      def self.querystring(classname="Endpoints")
        return CGI.escape("PREFIX ot:<http://www.opentox.org/api/1.1#>
        PREFIX ota:<http://www.opentox.org/algorithms.owl#>
        PREFIX owl:<http://www.w3.org/2002/07/owl#>
        PREFIX dc:<http://purl.org/dc/elements/1.1/>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX otee:<http://www.opentox.org/echaEndpoints.owl#>
        PREFIX toxcast:<http://www.opentox.org/toxcast.owl#>
        select *
          where {
            ?endpoint  rdfs:subClassOf  otee:#{classname}.
            ?endpoint dc:title ?title.
          }")
      end

      def self.make_option_list(endpoint="Endpoints", level=1)
      out = ""
        results = echa_endpoints(endpoint) rescue results = []
        results.each do |result|
          r = result.split(',')
          endpointname = r.first.split("#").last
          title = r[1..r.size-1]
          out += "<option value='#{r.first}' id='#{endpointname}' class='level_#{level}'>#{title}</option>\n"
          out += make_option_list(endpointname, level + 1)
        end
        return out
      end

      def self.endpoint_option_list(include_blank=true)
        out = "<select id='endpoint' name='endpoint'>\n"
        out += "<option value='' id='please_select'>Please select</option>\n" if include_blank 
        out += make_option_list
        out += "</select>\n"
        return out
      end

      # Gets Endpoints of specific level from ontology service
      # Top level with endpoint="Endpoints"
      # e.G. Ecotoxic effects endpoints with  endpoint="EcotoxicEffects"
      # if ontology service is not reachable it returns an empty array
      # @param [String] endpoint
      # @return [Array] of endpoints: e.G. "http://www.opentox.org/echaEndpoints.owl#EcotoxicEffects,Ecotoxic effects"
      def self.echa_endpoints(endpoint)
        begin
          RestClientWrapper.get("#{ONTOLOGY_SERVER}?query=#{querystring(endpoint)}",:accept => "text/csv").collect{|l| l.gsub("\r\n", "") if l.match(/^http/)}.uniq.compact.sort
        rescue
          LOGGER.warn "OpenTox::Ontology::Echa.echa_endpoints ontology service is not reachable."
          []
        end
      end

      def self.endpoints
        RestClientWrapper.get("http://apps.ideaconsult.net:8080/ambit2/query/ndatasets_endpoint",:accept => "text/csv").collect { |line| line.split(',').first if line.match(/^http/) }.compact
      end

      def self.datasets(endpoint)
        RestClientWrapper.get("http://apps.ideaconsult.net:8080/ambit2/dataset?feature_sameas=#{URI.encode endpoint}", :accept => "text/uri-list").split("\n")
      end

    end

    #Model Class for OpenTox::Ontology to register/deregister and check models in the ontology service
    #@example Register a model URI to the ontology service, check if model URI exists and delete it 
    #  uri = "http://mymodelservice.net/model/1" # model uri will be checked by the ontology service itself
    #  OpenTox::Ontology::Model.register(uri)
    #  puts OpenTox::Ontology::Model.exists?(uri) # => should return true
    #  OpenTox::Ontology::Model.delete(uri)
    #  puts OpenTox::Ontology::Model.exists?(uri) # => should return false
    module Model

      # Register an OpenTox resource into ontology service
      # @param [String] uri, URI of recource to register
      # @param [String] subjectid
      def self.register(uri, subjectid=nil)
        begin
          RestClientWrapper.post(ONTOLOGY_SERVER, {:uri => uri}, {:subjectid => CGI.escape(subjectid)})
        rescue
          LOGGER.warn "OpenTox::Ontology::Model.register ontology service is not reachable. Failed to register URI: #{uri} with subjectid: #{subjectid}"
          false
        end
      end

      # Deregister an OpenTox resource into ontology service
      # @param [String] uri, URI of recource to deregister/delete
      # @param [String] subjectid
      def self.delete(uri, subjectid=nil)
        begin
          RestClientWrapper.delete("#{ONTOLOGY_SERVER}?uri=#{CGI.escape(uri)}", {:subjectid => CGI.escape(subjectid)})
        rescue
          LOGGER.warn "OpenTox::Ontology::Model.exists ontology service is not reachable. Failed to delete URI: #{uri} with subjectid: #{subjectid}"
          false
        end
      end

      # Find an OpenTox resource in the ontology service
      # @param [String] uri, URI of recource to find
      # @param [String] subjectid
      def self.exists?(uri, subjectid=nil)
        begin
          out = RestClientWrapper.get("#{ONTOLOGY_SERVER}?query=#{querystring(uri)}",:accept => "text/csv").collect{|l| l.gsub("\r\n", "") if l.match(/^http/)}.uniq.compact
          return true if out.size > 0
          false
        rescue
          LOGGER.warn "OpenTox::Ontology::Model.exists ontology service is not reachable. Failed to check for URI: #{uri} with subjectid: #{subjectid}"
          false
        end
      end

      private
      # Query string to find a registered model
      # @param [String] uri, model URI
      def self.querystring(uri)
        return CGI.escape("PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX owl:<http://www.w3.org/2002/07/owl#>
        PREFIX dc:<http://purl.org/dc/elements/1.1/>
        PREFIX dcterms:<http://purl.org/dc/terms/>
        PREFIX ot:<http://www.opentox.org/api/1.1#>
        select distinct ?model ?title ?creator ?trainingDataset ?algorithm
          where {
          ?model rdf:type ot:Model;
            OPTIONAL {?model dc:title ?title}.
            OPTIONAL {?model dc:creator ?creator}.
            OPTIONAL {?model ot:trainingDataset ?trainingDataset}.
            OPTIONAL {?model ot:algorithm ?algorithm }.
            FILTER (?model = <#{uri}>)
          }")
      end
    end

  end
end