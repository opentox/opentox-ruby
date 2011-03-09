module OpenTox
  module OntologyService
    module Endpoints
      require 'sparql/client'
      @sparql = SPARQL::Client.new("http://apps.ideaconsult.net:8080/ontology")
      def self.qs(classname="Endpoints")
        return "PREFIX ot:<http://www.opentox.org/api/1.1#>
        PREFIX ota:<http://www.opentox.org/algorithms.owl#>
        PREFIX owl:<http://www.w3.org/2002/07/owl#>
        PREFIX dc:<http://purl.org/dc/elements/1.1/>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX otee:<http://www.opentox.org/echaEndpoints.owl#>
        PREFIX toxcast:<http://www.opentox.org/toxcast.owl#>
          select ?Endpoints ?title ?id
          where {?Endpoints rdfs:subClassOf otee:#{classname}.
          OPTIONAL {?Endpoints dc:title ?title}.
          OPTIONAL {?Endpoints dc:identifier ?id}.}
          ORDER BY ?title"
      end
      
      def self.make_option_list(endpoint="Endpoints", level=1)
        out = ""
        results = @sparql.query(qs(endpoint)) rescue results = []
        results.each do |result|
          endpointname = result.Endpoints.to_s.split('#').last
          title = result.bound?(:title) ? result.title : endpointname     
          out += "<option value='#{title}' id='#{endpointname}' class='level_#{level}'>#{title}</option>\n"
          out += make_option_list(endpointname, level + 1)
        end
        return out
      end
    
      def self.get_endpoint_selectlist(include_blank=true)
        out = "<select id='endpoint' name='endpoint'>\n"
        out += "<option value='' id='please_select'>Please select</option>\n" if include_blank 
        out += make_option_list
        out += "</select>\n"
        return out
      end
    end
  end
end