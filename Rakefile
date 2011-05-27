require 'rubygems'
require 'rake'

begin
  require 'jeweler'
  Jeweler::Tasks.new do |gem|
    gem.name = "opentox-ruby"
    gem.summary = %Q{Ruby wrapper for the OpenTox REST API}
    gem.description = %Q{Ruby wrapper for the OpenTox REST API (http://www.opentox.org)}
    gem.email = "helma@in-silico.ch"
    gem.homepage = "http://github.com/helma/opentox-ruby"
    gem.authors = ["Christoph Helma, Martin Guetlein, Andreas Maunz, Micha Rautenberg, David Vorgrimmler"]
    # dependencies with versions
    gem.add_dependency "sinatra", "=1.2.6"
    gem.add_dependency "emk-sinatra-url-for", "=0.2.1"
    gem.add_dependency "sinatra-respond_to", "=0.7.0"
    gem.add_dependency "sinatra-static-assets", "=0.5.0"
    gem.add_dependency "rest-client", "=1.6.1"
    gem.add_dependency "rack", "=1.3.0"
    gem.add_dependency "rack-contrib", "=1.1.0"
    gem.add_dependency "rack-flash", "=0.1.1"
    gem.add_dependency "nokogiri", "=1.4.4"
    gem.add_dependency "rubyzip", "=0.9.4"
    gem.add_dependency "roo", "=1.9.3"
    gem.add_dependency "spreadsheet", "=0.6.5.4"
    gem.add_dependency "google-spreadsheet-ruby", "=0.1.5"
    gem.add_dependency "yajl-ruby", "=0.8.2"
    gem.add_dependency "tmail", "=1.2.7.1"
    gem.add_dependency "rinruby", "=2.0.2"
    gem.add_dependency "ohm", "=0.1.3"
    gem.add_dependency "ohm-contrib", "=0.1.1"
    gem.add_dependency "SystemTimer", "=1.2.3"
    gem.add_dependency "rjb", "=1.3.4"
    gem.add_dependency "haml", "=3.1.1"
    #valiation-gems
    gem.add_dependency "dm-core",  "=1.1.0"
    gem.add_dependency "dm-serializer",  "=1.1.0"
    gem.add_dependency "dm-timestamps", "=1.1.0"
    gem.add_dependency "dm-types",  "=1.1.0"
    gem.add_dependency "dm-migrations",  "=1.1.0"
    gem.add_dependency "dm-validations",  "=1.1.0"
    gem.add_dependency "dm-sqlite-adapter", "=1.1.0"
    gem.add_dependency "ruby-plot", "=0.5.0"

    gem.add_development_dependency 'jeweler'
    gem.files =  FileList["[A-Z]*", "{bin,generators,lib,test}/**/*", 'lib/jeweler/templates/.gitignore']
  end
  Jeweler::GemcutterTasks.new
rescue LoadError
  puts "Jeweler (or a dependency) not available. Install it with: sudo gem install jeweler"
end

require 'rake/testtask'
Rake::TestTask.new(:test) do |test|
  test.libs << 'lib' << 'test'
  test.pattern = 'test/**/*_test.rb'
  test.verbose = true
end

begin
  require 'rcov/rcovtask'
  Rcov::RcovTask.new do |test|
    test.libs << 'test'
    test.pattern = 'test/**/*_test.rb'
    test.verbose = true
  end
rescue LoadError
  task :rcov do
    abort "RCov is not available. In order to run rcov, you must: sudo gem install spicycode-rcov"
  end
end

task :test => :check_dependencies

task :default => :test

require 'rake/rdoctask'
Rake::RDocTask.new do |rdoc|
  if File.exist?('VERSION')
    version = File.read('VERSION')
  else
    version = ""
  end

  rdoc.rdoc_dir = 'rdoc'
  rdoc.title = "opentox-ruby #{version}"
  rdoc.rdoc_files.include('README*')
  rdoc.rdoc_files.include('lib/**/*.rb')
end
