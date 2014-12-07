# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 'bindeps/version'

Gem::Specification.new do |spec|
  spec.name          = "bindeps"
  spec.version       = Bindeps::VERSION
  spec.authors       = ["Richard Smith-Unna", "Chris Boursnell"]
  spec.email         = ["rds45@cam.ac.uk"]
  spec.description   = %q{binary dependency management for ruby gems}
  spec.summary       = %q{binary dependency management for ruby gems}
  spec.homepage      = "https://github.com/Blahah/bindeps"
  spec.license       = "MIT"

  spec.files         = `git ls-files`.split($/)
  spec.executables   = spec.files.grep(%r{^bin/}) { |f| File.basename(f) }
  spec.test_files    = spec.files.grep(%r{^(test|spec|features)/})
  spec.require_paths = ["lib"]

  spec.add_dependency 'fixwhich', '~> 1.0', '>= 1.0.2'

  spec.add_development_dependency "bundler", "~> 1.3"
  spec.add_development_dependency 'rake', '~> 10.3', '>= 10.3.2'
  spec.add_development_dependency 'turn', '~> 0.9', '>= 0.9.7'
  spec.add_development_dependency 'minitest', '~> 4', '>= 4.7.5'
  spec.add_development_dependency 'simplecov', '~> 0.8', '>= 0.8.2'
  spec.add_development_dependency 'shoulda-context', '~> 1.2', '>= 1.2.1'
  spec.add_development_dependency 'coveralls', '~> 0.7'
end
