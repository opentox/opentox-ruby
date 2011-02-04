
# adding additional fields to Exception class to format errors according to OT-API
class Exception
  attr_accessor :errorCause
  def http_code; 500; end
end

module OpenTox
  
  class BadRequestError < RuntimeError
    def http_code; 400; end
  end
  
  class NotAuthorizedError < RuntimeError
    def http_code; 401; end
  end
  
  class NotFoundError < RuntimeError
    def http_code; 404; end
  end

  class ServiceUnavailableError < RuntimeError
    def http_code; 503; end
  end
  
  class RestCallError < RuntimeError
    attr_accessor :rest_params
    def http_code; 502; end
  end

  class ErrorReport
    
    # TODO replace params with URIs (errorCause -> OT.errorCause)
    attr_reader :message, :actor, :errorCause, :http_code, :errorDetails, :errorType

    private
    def initialize( http_code, erroType, message, actor, errorCause, rest_params=nil, backtrace=nil )
      @http_code = http_code
      @errorType = erroType
      @message = message
      @actor = actor
      @errorCause = errorCause
      @rest_params = rest_params
      @backtrace = backtrace
    end
    
    public
    # creates a error report object, from an ruby-exception object
    # @param [Exception] error
    # @param [String] actor, URI of the call that cause the error
    def self.create( error, actor )
      rest_params = error.rest_params if error.is_a?(OpenTox::RestCallError) and error.rest_params
      backtrace = error.backtrace.short_backtrace if CONFIG[:backtrace]
      ErrorReport.new( error.http_code, error.class.to_s, error.message, actor, error.errorCause, rest_params, backtrace )
    end
    
    def self.from_rdf(rdf)
      metadata = OpenTox::Parser::Owl.metadata_from_rdf( rdf, OT.ErrorReport )
      ErrorReport.new(metadata[OT.statusCode], metadata[OT.errorCode], metadata[OT.message], metadata[OT.actor], metadata[OT.errorCause])
    end
    
    # overwrite sorting to make easier readable
    def to_yaml_properties
       p = super
       p = ( p - ["@backtrace"]) + ["@backtrace"] if @backtrace
       p = ( p - ["@errorCause"]) + ["@errorCause"] if @errorCause
       p
    end
    
    def rdf_content()
      c = {
        RDF.type => OT.ErrorReport,
        OT.statusCode => @http_code,
        OT.message => @message,
        OT.actor => @actor,
        OT.errorCode => @errorType,
      }
      c[OT.errorCause] = @errorCause.rdf_content if @errorCause
      c
    end
    
    def to_rdfxml
      s = Serializer::Owl.new
      s.add_resource(CONFIG[:services]["opentox-task"]+"/tmpId/ErrorReport/tmpId", OT.errorReport, rdf_content)
      s.to_rdfxml
    end
  end
end

class Array
  def short_backtrace
    short = []
    each do |c|
      break if c =~ /sinatra\/base/
      short << c
    end
    short.join("\n")
  end
end