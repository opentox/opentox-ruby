
OT_LOGO = "http://opentox.informatik.uni-freiburg.de/ot-logo.png"


class String
  
  # encloses URI in text with with link tag
  # @return [String] new text with marked links
  def link_urls
    self.gsub(/(?i)http(s?):\/\/[^\r\n\s']*/, '<a href=\0>\0</a>')
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
  def self.text_to_html( text, related_links=nil, description=nil, post_params=nil )
    
    # TODO add title as parameter
    title = nil #$sinatra.url_for($sinatra.request.env['PATH_INFO'], :full) if $sinatra
    
    html = <<EOF
<html>
EOF
    html.chomp!
    html += "<title>"+title+"</title>" if title
    html += <<EOF 
<img src="
EOF
    html.chomp!
    html += OT_LOGO
    html += <<EOF 
">
<body>
EOF
   html.chomp!
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
   html += <<EOF
<pre>
<p style="padding:15px; border:10px solid #5D308A">
EOF
   html.chomp!
   html += text.link_urls
   html += <<EOF
</p>
</pre>
</body>
<html>
EOF
    html
  end
  
end

#puts OpenTox.text_to_html("bla")