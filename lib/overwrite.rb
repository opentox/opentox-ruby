# class overwrites aka monkey patches
# hack: store sinatra instance in global var $url_provider to make url_for and halt methods accessible
before {
  raise "should not happen, url provider already differently initialized "+
    $url_provider.request.host.to_s+" != "+self.request.host.to_s if
    $url_provider and $url_provider.request.host!=self.request.host and 
    $url_provider.request.script_name!=self.request.script_name
  $url_provider = self
  # stupid internet explorer does not ask for text/html, add this manually 
  request.env['HTTP_ACCEPT'] += ";text/html" if request.env["HTTP_USER_AGENT"]=~/MSIE/
  request.env['HTTP_ACCEPT']=request.params["media"] if request.params["media"]
}


class Ohm::Model
  def self.find(p)
    if (p!=nil)
      p.delete(:captures)
      p.delete(:splat)
    end
    super p
  end
end   

  
# Error handling
# Errors are logged as error and formated according to acccept-header
# Non OpenTox::Errors (defined in error.rb) are handled as internal error (500), stacktrace is logged
# IMPT: set sinatra settings :show_exceptions + :raise_errors to false in config.ru, otherwise Rack::Showexceptions takes over
error Exception do
  error = request.env['sinatra.error']
  # log error message and backtrace to logfile
  LOGGER.error error.class.to_s+": "+error.message
  LOGGER.error ":\n"+error.backtrace.join("\n")
  
  actor = "#{request.env['rack.url_scheme']}://#{request.env['HTTP_HOST']}#{request.env['REQUEST_URI']}"
  rep = OpenTox::ErrorReport.create(error, actor)

  case request.env['HTTP_ACCEPT']
  when /rdf/
    content_type 'application/rdf+xml'
    halt error.http_code,rep.to_rdfxml
  when /html/
    content_type 'text/html'
    halt error.http_code,(OpenTox.text_to_html rep.to_yaml, @subjectid)
  else
    content_type 'application/x-yaml'
    halt error.http_code,rep.to_yaml
  end
end

class Sinatra::Base

  def return_task( task )
    raise "http_code == nil" unless task.http_code!=nil
    case request.env['HTTP_ACCEPT']
    when /rdf/
      response['Content-Type'] = "application/rdf+xml"
      halt task.http_code,task.to_rdfxml
    when /yaml/
      response['Content-Type'] = "application/x-yaml"
      halt task.http_code,task.to_yaml # PENDING differs from task-webservice
    when /html/
      response['Content-Type'] = "text/html"
      # html -> task created with html form -> redirect to task uri
      redirect task.uri 
    else # default /uri-list/
      response['Content-Type'] = "text/uri-list"
      if task.completed?
        halt task.http_code,task.resultURI+"\n"
      else
        halt task.http_code,task.uri+"\n"
      end
    end
  end
end

class String
  def task_uri?
    self.uri? && !self.match(/task/).nil?
  end
  
  def dataset_uri?
   self.uri? && !self.match(/dataset/).nil?
  end
 
  def self.model_uri?
   self.uri? && !self.match(/model/).nil?
  end

  def uri?
    begin
      u = URI::parse(self)
      return (u.scheme!=nil and u.host!=nil)
    rescue URI::InvalidURIError
      return false
    end
  end

  def underscore
    self.gsub(/::/, '/').
    gsub(/([A-Z]+)([A-Z][a-z])/,'\1_\2').
    gsub(/([a-z\d])([A-Z])/,'\1_\2').
    tr("-", "_").
    downcase
  end
end

require 'logger'
# logging
#class Logger
class OTLogger < Logger
  
  def pwd
    path = Dir.pwd.to_s
    index = path.rindex(/\//)
    return path if index==nil
    path[(index+1)..-1]
  end
  
  def trace()
    lines = caller(0)
    n = 2
    line = lines[n]
    
    while (line =~ /spork.rb/ or line =~ /create/ or line =~ /overwrite.rb/)
      n += 1
      line = lines[n]
    end
  
    index = line.rindex(/\/.*\.rb/)
    return line if index==nil
    line[index..-1]
  end
  
  def format(msg)
    pwd.ljust(18)+" :: "+msg.to_s+"           :: "+trace
  end
  
  def debug(msg)
    super format(msg)
  end
  
  def info(msg)
    super format(msg)
  end
  
  def warn(msg)
    super format(msg)
  end

  def error(msg)
    super format(msg)
  end

end

# make migration from datamapper more straightforward
class Ohm::Model
  def self.get(id)
    self[id]
  end
end
