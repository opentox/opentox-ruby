# Generated by jeweler
# DO NOT EDIT THIS FILE DIRECTLY
# Instead, edit Jeweler::Tasks in Rakefile, and run 'rake gemspec'
# -*- encoding: utf-8 -*-

Gem::Specification.new do |s|
  s.name = %q{opentox-ruby}
  s.version = "3.1.0"

  s.required_rubygems_version = Gem::Requirement.new(">= 0") if s.respond_to? :required_rubygems_version=
  s.authors = ["Christoph Helma, Martin Guetlein, Andreas Maunz, Micha Rautenberg, David Vorgrimmler"]
  s.date = %q{2012-03-13}
  s.description = %q{Ruby wrapper for the OpenTox REST API (http://www.opentox.org)}
  s.email = %q{helma@in-silico.ch}
  s.extra_rdoc_files = [
    "ChangeLog",
    "LICENSE",
    "README.markdown"
  ]
  s.files = [
    "ChangeLog",
    "LICENSE",
    "README.markdown",
    "Rakefile",
    "VERSION",
    "lib/algorithm.rb",
    "lib/authorization.rb",
    "lib/compound.rb",
    "lib/config/config_ru.rb",
    "lib/dataset.rb",
    "lib/environment.rb",
    "lib/error.rb",
    "lib/feature.rb",
    "lib/helper.rb",
    "lib/model.rb",
    "lib/ontology.rb",
    "lib/opentox-ruby.rb",
    "lib/opentox.owl",
    "lib/opentox.rb",
    "lib/overwrite.rb",
    "lib/parser.rb",
    "lib/policy.rb",
    "lib/r-util.rb",
    "lib/rest_client_wrapper.rb",
    "lib/serializer.rb",
    "lib/spork.rb",
    "lib/stratification.R",
    "lib/task.rb",
    "lib/templates/default_guest_policy.xml",
    "lib/templates/default_policy.xml",
    "lib/to-html.rb",
    "lib/transform.rb",
    "lib/utils.rb",
    "lib/validation.rb"
  ]
  s.homepage = %q{http://github.com/opentox/opentox-ruby}
  s.require_paths = ["lib"]
  s.rubygems_version = %q{1.5.3}
  s.summary = %q{Ruby wrapper for the OpenTox REST API}

  if s.respond_to? :specification_version then
    s.specification_version = 3

    if Gem::Version.new(Gem::VERSION) >= Gem::Version.new('1.2.0') then
      s.add_runtime_dependency(%q<sinatra>, ["= 1.2.6"])
      s.add_runtime_dependency(%q<emk-sinatra-url-for>, ["= 0.2.1"])
      s.add_runtime_dependency(%q<sinatra-respond_to>, ["= 0.7.0"])
      s.add_runtime_dependency(%q<sinatra-static-assets>, ["= 0.5.0"])
      s.add_runtime_dependency(%q<rest-client>, ["= 1.6.1"])
      s.add_runtime_dependency(%q<rack>, ["= 1.3.5"])
      s.add_runtime_dependency(%q<rack-contrib>, ["= 1.1.0"])
      s.add_runtime_dependency(%q<rack-flash>, ["= 0.1.1"])
      s.add_runtime_dependency(%q<nokogiri>, ["= 1.4.4"])
      s.add_runtime_dependency(%q<rubyzip>, ["= 0.9.4"])
      s.add_runtime_dependency(%q<roo>, ["= 1.9.3"])
      s.add_runtime_dependency(%q<spreadsheet>, ["= 0.6.5.4"])
      s.add_runtime_dependency(%q<google-spreadsheet-ruby>, ["= 0.1.5"])
      s.add_runtime_dependency(%q<yajl-ruby>, ["= 0.8.2"])
      s.add_runtime_dependency(%q<rinruby>, ["= 2.0.2"])
      s.add_runtime_dependency(%q<ohm>, ["= 0.1.3"])
      s.add_runtime_dependency(%q<ohm-contrib>, ["= 0.1.1"])
      s.add_runtime_dependency(%q<SystemTimer>, ["= 1.2.3"])
      s.add_runtime_dependency(%q<rjb>, ["= 1.3.4"])
      s.add_runtime_dependency(%q<haml>, ["= 3.1.1"])
      s.add_runtime_dependency(%q<akephalos>, ["= 0.2.5"])
      s.add_runtime_dependency(%q<dm-core>, ["= 1.1.0"])
      s.add_runtime_dependency(%q<dm-serializer>, ["= 1.1.0"])
      s.add_runtime_dependency(%q<dm-timestamps>, ["= 1.1.0"])
      s.add_runtime_dependency(%q<dm-types>, ["= 1.1.0"])
      s.add_runtime_dependency(%q<dm-migrations>, ["= 1.1.0"])
      s.add_runtime_dependency(%q<dm-validations>, ["= 1.1.0"])
      s.add_runtime_dependency(%q<dm-sqlite-adapter>, ["= 1.1.0"])
      s.add_runtime_dependency(%q<ruby-plot>, ["= 0.6.0"])
      s.add_runtime_dependency(%q<gsl>, ["= 1.14.7"])
      s.add_runtime_dependency(%q<statsample>, ["= 1.1.0"])
      s.add_development_dependency(%q<jeweler>, [">= 0"])
    else
      s.add_dependency(%q<sinatra>, ["= 1.2.6"])
      s.add_dependency(%q<emk-sinatra-url-for>, ["= 0.2.1"])
      s.add_dependency(%q<sinatra-respond_to>, ["= 0.7.0"])
      s.add_dependency(%q<sinatra-static-assets>, ["= 0.5.0"])
      s.add_dependency(%q<rest-client>, ["= 1.6.1"])
      s.add_dependency(%q<rack>, ["= 1.3.5"])
      s.add_dependency(%q<rack-contrib>, ["= 1.1.0"])
      s.add_dependency(%q<rack-flash>, ["= 0.1.1"])
      s.add_dependency(%q<nokogiri>, ["= 1.4.4"])
      s.add_dependency(%q<rubyzip>, ["= 0.9.4"])
      s.add_dependency(%q<roo>, ["= 1.9.3"])
      s.add_dependency(%q<spreadsheet>, ["= 0.6.5.4"])
      s.add_dependency(%q<google-spreadsheet-ruby>, ["= 0.1.5"])
      s.add_dependency(%q<yajl-ruby>, ["= 0.8.2"])
      s.add_dependency(%q<rinruby>, ["= 2.0.2"])
      s.add_dependency(%q<ohm>, ["= 0.1.3"])
      s.add_dependency(%q<ohm-contrib>, ["= 0.1.1"])
      s.add_dependency(%q<SystemTimer>, ["= 1.2.3"])
      s.add_dependency(%q<rjb>, ["= 1.3.4"])
      s.add_dependency(%q<haml>, ["= 3.1.1"])
      s.add_dependency(%q<akephalos>, ["= 0.2.5"])
      s.add_dependency(%q<dm-core>, ["= 1.1.0"])
      s.add_dependency(%q<dm-serializer>, ["= 1.1.0"])
      s.add_dependency(%q<dm-timestamps>, ["= 1.1.0"])
      s.add_dependency(%q<dm-types>, ["= 1.1.0"])
      s.add_dependency(%q<dm-migrations>, ["= 1.1.0"])
      s.add_dependency(%q<dm-validations>, ["= 1.1.0"])
      s.add_dependency(%q<dm-sqlite-adapter>, ["= 1.1.0"])
      s.add_dependency(%q<ruby-plot>, ["= 0.6.0"])
      s.add_dependency(%q<gsl>, ["= 1.14.7"])
      s.add_dependency(%q<statsample>, ["= 1.1.0"])
      s.add_dependency(%q<jeweler>, [">= 0"])
    end
  else
    s.add_dependency(%q<sinatra>, ["= 1.2.6"])
    s.add_dependency(%q<emk-sinatra-url-for>, ["= 0.2.1"])
    s.add_dependency(%q<sinatra-respond_to>, ["= 0.7.0"])
    s.add_dependency(%q<sinatra-static-assets>, ["= 0.5.0"])
    s.add_dependency(%q<rest-client>, ["= 1.6.1"])
    s.add_dependency(%q<rack>, ["= 1.3.5"])
    s.add_dependency(%q<rack-contrib>, ["= 1.1.0"])
    s.add_dependency(%q<rack-flash>, ["= 0.1.1"])
    s.add_dependency(%q<nokogiri>, ["= 1.4.4"])
    s.add_dependency(%q<rubyzip>, ["= 0.9.4"])
    s.add_dependency(%q<roo>, ["= 1.9.3"])
    s.add_dependency(%q<spreadsheet>, ["= 0.6.5.4"])
    s.add_dependency(%q<google-spreadsheet-ruby>, ["= 0.1.5"])
    s.add_dependency(%q<yajl-ruby>, ["= 0.8.2"])
    s.add_dependency(%q<rinruby>, ["= 2.0.2"])
    s.add_dependency(%q<ohm>, ["= 0.1.3"])
    s.add_dependency(%q<ohm-contrib>, ["= 0.1.1"])
    s.add_dependency(%q<SystemTimer>, ["= 1.2.3"])
    s.add_dependency(%q<rjb>, ["= 1.3.4"])
    s.add_dependency(%q<haml>, ["= 3.1.1"])
    s.add_dependency(%q<akephalos>, ["= 0.2.5"])
    s.add_dependency(%q<dm-core>, ["= 1.1.0"])
    s.add_dependency(%q<dm-serializer>, ["= 1.1.0"])
    s.add_dependency(%q<dm-timestamps>, ["= 1.1.0"])
    s.add_dependency(%q<dm-types>, ["= 1.1.0"])
    s.add_dependency(%q<dm-migrations>, ["= 1.1.0"])
    s.add_dependency(%q<dm-validations>, ["= 1.1.0"])
    s.add_dependency(%q<dm-sqlite-adapter>, ["= 1.1.0"])
    s.add_dependency(%q<ruby-plot>, ["= 0.6.0"])
    s.add_dependency(%q<gsl>, ["= 1.14.7"])
    s.add_dependency(%q<statsample>, ["= 1.1.0"])
    s.add_dependency(%q<jeweler>, [">= 0"])
  end
end

