['rubygems', 'sinatra', 'sinatra/url_for', 'ohm', 'rest_client', 'yaml', 'cgi', 'spork', 'error', 'overwrite', 'environment'].each do |lib|
	require lib
end

begin
	require 'openbabel'
rescue LoadError
	puts "Please install Openbabel with 'rake openbabel:install' in the compound component"
end

['opentox', 'compound','dataset', 'parser','serializer', 'algorithm','model','task','validation','feature', 
 'rest_client_wrapper', 'authorization', 'policy', 'helper', 'to-html', 'ontology' ].each do |lib|
	require lib
end
