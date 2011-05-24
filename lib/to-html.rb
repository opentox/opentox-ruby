
OT_LOGO = File.join(CONFIG[:services]["opentox-validation"],"resources/ot-logo.png")

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
  # @example post params:
  # [ [ [:mandatory_param_1], [:mandatory_param_2], [:optional_param,"default_value"] ],
  #   [ [:alteranative_mandatory_param_1], [:alteranative_mandatory_param_2] ]
  # ] 
  # @param [String] text this is the actual content, 
  # @param [optional,String] related_links info on related resources
  # @param [optional,String] description general info
  # @param [optional,Array] post_params, array of arrays containing info on POST operation, see example
  # @return [String] html page
  def self.text_to_html( text, subjectid=nil, related_links=nil, description=nil, post_params=nil  )
    
    # TODO add title as parameter
    title = nil #$sinatra.url_for($sinatra.request.env['PATH_INFO'], :full) if $sinatra
    html = "<html>"
    html += "<title>"+title+"</title>" if title
    html += "<img src=\""+OT_LOGO+"\"><\/img><body>"
      
    if AA_SERVER
      user = OpenTox::Authorization.get_user(subjectid) if subjectid
      html +=  "<pre><p align=\"right\">"
      unless user
        html += "You are currently not logged in to "+$url_provider.url_for("",:full)+
          ", <a href="+$url_provider.url_for("/login",:full)+">login</a>"
      else
        html += "You are logged in as '#{user}' to "+$url_provider.url_for("",:full)+
          ", <a href="+$url_provider.url_for("/logout",:full)+">logout</a>"
      end
      html += "  </p></pre>"
    end 
   
    html += "<h3>Description</h3><pre><p>"+description.link_urls+"</p></pre>" if description
    html += "<h3>Related links</h3><pre><p>"+related_links.link_urls+"</p></pre>" if related_links
    if post_params
      html += "<h3>POST parameters</h3>"
      count = 0
      post_params.each do |p|
        html += "<pre><p>alternatively:</p></pre>" if count > 0
        html += "<pre><p><table><thead><tr><th>param</th><th>default_value</th></tr></thead>"
        p.each do |k,v|
          html += "<tr><th>"+k.to_s+"</th><th>"+(v!=nil ? v.to_s : "<i>mandatory</i>")+"</th></tr>"
        end
        html += "</table></p></pre>"
        count += 1
      end
    end
    html += "<h3>Content</h3>" if description || related_links
    html += "<pre><p style=\"padding:15px; border:10px solid \#5D308A\">"
    html += text.link_urls
    html += "</p></pre></body></html>"
    html
  end
  
  def self.login( msg=nil )
    html = "<html><title>Login</title><img src="+OT_LOGO+"><body>"
    html += "<form method='POST' action='"+$url_provider.url_for("/login",:full)+"'>"
    html += "<pre><p style=\"padding:15px; border:10px solid \#5D308A\">"
    html += msg+"\n\n" if msg
    html += "Please login to "+$url_provider.url_for("",:full)+"\n\n"
    html += "<table border=0>"
    html += "<tr><td>user:</td><td><input type='text' name='user' size='15' /></td></tr>"+
          "<tr><td>password:</td><td><input type='password' name='password' size='15' /></td></tr>"+
          #"<input type=hidden name=back_to value="+back_to.to_s+">"+
          "<tr><td><input type='submit' value='Login' /></td></tr>"
    html += "</table></p></pre></form></body></html>"
    html
  end
end

=begin
get '/logout/?' do
  response.set_cookie("subjectid",{:value=>nil})
  content_type "text/html"
  content = "Sucessfully logged out from "+$url_provider.url_for("",:full)
  OpenTox.text_to_html(content)
end

get '/login/?' do
  content_type "text/html"
  OpenTox.login
end

post '/login/?' do
  subjectid = OpenTox::Authorization.authenticate(params[:user], params[:password])
  if (subjectid)
    response.set_cookie("subjectid",{:value=>subjectid})
    content_type "text/html"
    content = "Sucessfully logged in as '"+params[:user]+"' to "+$url_provider.url_for("",:full)
    OpenTox.text_to_html(content,subjectid)    
  else
    content_type "text/html"
    OpenTox.login("Login failed, please try again")
  end
end
=end

