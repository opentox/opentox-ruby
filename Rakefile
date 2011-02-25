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
    # dependencies
    [ "sinatra",
      "emk-sinatra-url-for",
      "sinatra-respond_to",
      "sinatra-static-assets",
      "rest-client",
      "rack",
      "rack-contrib",
      "rack-flash",
      "nokogiri",
      "rubyzip",
      "roo",
      "spreadsheet",
      "google-spreadsheet-ruby",
      "yajl-ruby",
      "tmail",
      "rinruby",
      "ohm",
      "rjb"
    ].each { |dep| gem.add_dependency dep }
    [ "dm-core",
      'dm-serializer',
      'dm-timestamps',
      'dm-types',
      'dm-migrations',
      "dm-mysql-adapter",
      "dm-validations",
    ].each {|dep| gem.add_dependency dep, ">= 1" }
    gem.add_dependency "haml", ">=3"
    ['jeweler'].each { |dep| gem.add_development_dependency dep }
    gem.files =  FileList["[A-Z]*", "{bin,generators,lib,test}/**/*", 'lib/jeweler/templates/.gitignore']
    #gem.files.include %w(lib/environment.rb, lib/algorithm.rb, lib/compound.rb, lib/dataset.rb, lib/model.rb, lib/validation.rb, lib/templates/*)
    # gem is a Gem::Specification... see http://www.rubygems.org/read/chapter/20 for additional settings
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
