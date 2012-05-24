
OT_LOGO = "/" + CONFIG[:services]["opentox-validation"].split("/").last + "/resources/ot-logo.png"

class String
  
  # encloses URI in text with with link tag
  # @return [String] new text with marked links
  def link_urls
    self.gsub(/(?i)http(s?):\/\/[^\r\n\s']*/, '<a href="\0">\0</a>')
  end
end

module OpenTox
  
  # produces a html page for making web services browser friendly
  # format of text (=string params) is preserved (e.g. line breaks)
  # urls are marked as links
  #
  # @param [String] text this is the actual content, 
  # @param [optional,String] related_links info on related resources
  # @param [optional,String] description general info
  # @param [optional,Array] post_command, infos for the post operation, object defined below
  # @return [String] html page
  def self.text_to_html( text, subjectid=nil, related_links=nil, description=nil, post_command=nil  )
    
    # TODO add title as parameter
    title = nil #$sinatra.url_for($sinatra.request.env['PATH_INFO'], :full) if $sinatra
    html = "<html>"
    html += "<title>"+title+"</title>" if title
    html += "<img src=\""+OT_LOGO+"\"><\/img><body>"
      
    if AA_SERVER
      user = OpenTox::Authorization.get_user(subjectid) if subjectid
      html +=  "<pre><p align=\"right\">"
      unless user
        html += "You are currently not signed in to "+$url_provider.request.host.to_s+
          ", <a href="+$url_provider.url_for("/sign_in",:full)+">sign in</a>"
      else
        html += "You are signed in as '#{user}' to "+$url_provider.request.host.to_s+
          ", <a href="+$url_provider.url_for("/sign_out",:full)+">sign out</a>"
      end
      html += "  </p></pre>"
    end 
   
    html += "<h3>Description</h3><pre><p>"+description.link_urls+"</p></pre>" if description
    html += "<h3>Related links</h3><pre><p>"+related_links.link_urls+"</p></pre>" if related_links
    if post_command
      raise "not a post command" unless post_command.is_a?(OpenTox::PostCommand)
      html += "<h3>POST command</h3>"
      html += post_command.to_html
    end
    html += "<h3>Content</h3>" if description || related_links || post_command
    html += "<pre><p style=\"padding:15px; border:10px solid \#5D308A\">"
    html += text.link_urls
    html += "</p></pre></body></html>"
    html
  end
  
  def self.sign_in( msg=nil )
    html = "<html><title>Login</title><img src="+OT_LOGO+"><body>"
    html += "<form method='POST' action='"+$url_provider.url_for("/sign_in",:full)+"'>"
    html += "<pre><p style=\"padding:15px; border:10px solid \#5D308A\">"
    html += msg+"\n\n" if msg
    html += "Please sign in to "+$url_provider.request.host.to_s+"\n\n"
    html += "<table border=0>"
    html += "<tr><td>user:</td><td><input type='text' name='user' size='15' /></td></tr>"+
          "<tr><td>password:</td><td><input type='password' name='password' size='15' /></td></tr>"+
          #"<input type=hidden name=back_to value="+back_to.to_s+">"+
          "<tr><td><input type='submit' value='Sign in' /></td></tr>"
    html += "</table></p></pre></form></body></html>"
    html
  end
  
  class PostAttribute
    attr_accessor :name, :is_mandatory, :default, :description
    
    def initialize(name, is_mandatory=true, default=nil, description=nil)
      @name = name
      @is_mandatory = is_mandatory
      @default = default
      @description = description
    end
  end
  
  class PostCommand
    attr_accessor :attributes, :uri, :name
    
    def initialize( uri, name="Send" )
      @uri = uri
      @name = name
      @attributes = []
    end
   
    def to_html
      html = "<form method='POST' action='"+@uri.to_s+"'>"
      html << "<pre><p>"
      html << "<table border=0>"
      #html << "<tr><td colspan='3'><i><sup>Mandatory params are marked with *.</sup></i></td></tr>"
      attributes.each do |a|
        mandatory_string = a.is_mandatory ? "*" : ""
        html << "<tr><td>"+a.name.to_s+":"+mandatory_string+"</td>"
        html << "<td><input type='text' name='"+a.name.to_s+
          "' size='50' value='"+a.default.to_s+"'/></td>"
        html << "<td><i><sup>"+a.description.to_s+"</sup></i></td></tr>"
      end
      html << "<tr><td colspan='3'><input type='submit' value='"+@name.to_s+"' /></td></tr>"
      html << "</table></p></pre></form>"
      html
    end
  end
end

get '/sign_out/?' do
  logout
  content_type "text/html"
  content = "Sucessfully signed out from "+$url_provider.request.host.to_s+" ( Back to "+
      $url_provider.url_for("",:full)+" )"
  OpenTox.text_to_html(content)
end

get '/sign_in/?' do
  content_type "text/html"
  OpenTox.sign_in
end

post '/sign_in/?' do
  subjectid = login(params[:user], params[:password])
  if (subjectid)
    content_type "text/html"
    content = "Sucessfully signed in as '"+params[:user]+"' to "+$url_provider.request.host.to_s+" ( Back to "+
      $url_provider.url_for("",:full)+" )"
    OpenTox.text_to_html(content,subjectid)    
  else
    content_type "text/html"
    OpenTox.sign_in("Login failed, please try again")
  end
end

