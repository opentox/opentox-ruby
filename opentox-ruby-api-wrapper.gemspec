# Generated by jeweler
# DO NOT EDIT THIS FILE
# Instead, edit Jeweler::Tasks in Rakefile, and run `rake gemspec`
# -*- encoding: utf-8 -*-

Gem::Specification.new do |s|
  s.name = %q{opentox-ruby-api-wrapper}
  s.version = "1.2.3"

  s.required_rubygems_version = Gem::Requirement.new(">= 0") if s.respond_to? :required_rubygems_version=
  s.authors = ["Christoph Helma"]
  s.date = %q{2009-12-21}
  s.description = %q{Ruby wrapper for the OpenTox REST API (http://www.opentox.org)}
  s.email = %q{helma@in-silico.ch}
  s.executables = ["opentox-install-debian.sh", "yaml2owl.rb"]
  s.extra_rdoc_files = [
    "LICENSE",
     "README.rdoc"
  ]
  s.files = [
    "LICENSE",
     "README.rdoc",
     "Rakefile",
     "VERSION",
     "bin/opentox-install-debian.sh",
     "bin/yaml2owl.rb",
     "lib/algorithm.rb",
     "lib/compound.rb",
     "lib/dataset.rb",
     "lib/environment.rb",
     "lib/helper.rb",
     "lib/model.rb",
     "lib/opentox-ruby-api-wrapper.rb",
     "lib/opentox.owl",
     "lib/owl.rb",
     "lib/spork.rb",
     "lib/task.rb",
     "lib/tasks/opentox.rb",
     "lib/tasks/redis.rb",
     "lib/templates/config.ru",
     "lib/templates/config.ru",
     "lib/templates/config.yaml",
     "lib/templates/config.yaml",
     "lib/utils.rb"
  ]
  s.has_rdoc = true
  s.homepage = %q{http://github.com/helma/opentox-ruby-api-wrapper}
  s.rdoc_options = ["--charset=UTF-8"]
  s.require_paths = ["lib"]
  s.rubygems_version = %q{1.3.1}
  s.summary = %q{Ruby wrapper for the OpenTox REST API}

  if s.respond_to? :specification_version then
    current_version = Gem::Specification::CURRENT_SPECIFICATION_VERSION
    s.specification_version = 2

    if Gem::Version.new(Gem::RubyGemsVersion) >= Gem::Version.new('1.2.0') then
      s.add_runtime_dependency(%q<rest-client>, [">= 0"])
      s.add_runtime_dependency(%q<sinatra>, [">= 0"])
      s.add_runtime_dependency(%q<rack>, [">= 0"])
      s.add_runtime_dependency(%q<rack-contrib>, [">= 0"])
      s.add_runtime_dependency(%q<thin>, [">= 0"])
      s.add_runtime_dependency(%q<emk-sinatra-url-for>, [">= 0"])
      s.add_runtime_dependency(%q<cehoffman-sinatra-respond_to>, [">= 0"])
      s.add_development_dependency(%q<cucumber>, [">= 0"])
    else
      s.add_dependency(%q<rest-client>, [">= 0"])
      s.add_dependency(%q<sinatra>, [">= 0"])
      s.add_dependency(%q<rack>, [">= 0"])
      s.add_dependency(%q<rack-contrib>, [">= 0"])
      s.add_dependency(%q<thin>, [">= 0"])
      s.add_dependency(%q<emk-sinatra-url-for>, [">= 0"])
      s.add_dependency(%q<cehoffman-sinatra-respond_to>, [">= 0"])
      s.add_dependency(%q<cucumber>, [">= 0"])
    end
  else
    s.add_dependency(%q<rest-client>, [">= 0"])
    s.add_dependency(%q<sinatra>, [">= 0"])
    s.add_dependency(%q<rack>, [">= 0"])
    s.add_dependency(%q<rack-contrib>, [">= 0"])
    s.add_dependency(%q<thin>, [">= 0"])
    s.add_dependency(%q<emk-sinatra-url-for>, [">= 0"])
    s.add_dependency(%q<cehoffman-sinatra-respond_to>, [">= 0"])
    s.add_dependency(%q<cucumber>, [">= 0"])
  end
end
