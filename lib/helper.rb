helpers do

  # Authentification
  def protected!(subjectid)
    if env["session"]
      unless authorized?(subjectid)
        flash[:notice] = "You don't have access to this section: "
        redirect back
      end
    elsif !env["session"] && subjectid
      unless authorized?(subjectid)
        LOGGER.debug "URI not authorized: clean: " + clean_uri("#{request.env['rack.url_scheme']}://#{request.env['HTTP_HOST']}#{request.env['REQUEST_URI']}").to_s + " full: #{request.env['rack.url_scheme']}://#{request.env['HTTP_HOST']}#{request.env['REQUEST_URI']} with request: #{request.env['REQUEST_METHOD']}"
        raise OpenTox::NotAuthorizedError.new "Not authorized" 
      end
    else
      raise OpenTox::NotAuthorizedError.new "Not authorized" unless authorized?(subjectid)
    end
  end

  #Check Authorization for URI with method and subjectid. 
  def authorized?(subjectid)
    # hack for reports, address problem as soon as subjectid is not longer allowed as param 
    return true if request.env['REQUEST_URI'] =~ /validation\/report\/.*svg$/
    request_method = request.env['REQUEST_METHOD']
    uri = clean_uri("#{request.env['rack.url_scheme']}://#{request.env['HTTP_HOST']}#{request.env['REQUEST_URI']}")
    request_method = "GET" if request_method == "POST" &&  uri =~ /\/model\/\d+\/?$/
    return OpenTox::Authorization.authorized?(uri, request_method, subjectid)
  end

  #cleans URI from querystring and file-extension. Sets port 80 to emptystring
  # @param [String] uri 
  def clean_uri(uri)
    uri = uri.sub(" ", "%20")          #dirty hacks => to fix
    uri = uri[0,uri.index("InChI=")] if uri.index("InChI=") 
    
    out = URI.parse(uri)
    out.path = out.path[0, out.path.length - (out.path.reverse.rindex(/\/{1}\d+\/{1}/))] if out.path.index(/\/{1}\d+\/{1}/)  #cuts after /id/ for a&a 
    "#{out.scheme}:" + (out.port != 80 ? out.port : "") + "//#{out.host}#{out.path.chomp('/')}"
  end

  #unprotected uri for login
  def login_requests
    return env['REQUEST_URI'] =~ /\/login$/ 
   end

end

before do
  unless !AA_SERVER or login_requests or CONFIG[:authorization][:free_request].include?(env['REQUEST_METHOD']) 
    begin
      subjectid = nil
      subjectid = session[:subjectid] if session[:subjectid]
      subjectid = params[:subjectid]  if params[:subjectid] and !subjectid
      subjectid = request.env['HTTP_SUBJECTID'] if request.env['HTTP_SUBJECTID'] and !subjectid
      # see http://rack.rubyforge.org/doc/SPEC.html
      subjectid = CGI.unescape(subjectid) if subjectid.include?("%23")
      @subjectid = subjectid
    rescue
      #LOGGER.debug "OpenTox ruby api wrapper: helper before filter: NO subjectid for URI: #{request.env['rack.url_scheme']}://#{request.env['HTTP_HOST']}#{request.env['REQUEST_URI']}"
      subjectid = ""
    end
    @subjectid = subjectid
    protected!(subjectid)
  end
end

