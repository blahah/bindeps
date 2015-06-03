# Bindeps

[![Build Status](https://travis-ci.org/Blahah/bindeps.svg)](https://travis-ci.org/Blahah/bindeps)
[![Dependency Status](http://img.shields.io/gemnasium/Blahah/bindeps.svg)](https://gemnasium.com/Blahah/bindeps)
[![Coverage Status](http://img.shields.io/coveralls/Blahah/bindeps.svg)](https://coveralls.io/r/Blahah/bindeps)
[![Code Climate](http://img.shields.io/codeclimate/github/Blahah/bindeps.svg)](https://codeclimate.com/github/Blahah/bindeps)
[![Gem Version](http://img.shields.io/gem/v/bindeps.svg)](https://rubygems.org/gems/bindeps)
[![License](http://img.shields.io/:license-mit-blue.svg)](http://Blahah.mit-license.org)

Simple binary dependency management for Ruby gems

## Installation

```bash
$ gem install bindeps
```

## Usage

Create a YAML file describing your dependencies as a dictionary. Read the [bindeps YAML format specifications](https://github.com/Blahah/bindeps/wiki/bindeps_YAML_format_specifications).

```yaml
bowtie2:
  binaries:
    - bowtie2
    - bowtie2-align-l
    - bowtie2-align-s
    - bowtie2-build
    - bowtie2-build-l
    - bowtie2-build-s
    - bowtie2-inspect
    - bowtie2-inspect-l
    - bowtie2-inspect-s
  version:
    number: '2.2.3'
    command: 'bowtie2 --version'
  url:
    64bit:
      linux: http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.2.3/bowtie2-2.2.3-linux-x86_64.zip
      macosx: http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.2.3/bowtie2-2.2.3-macos-x86_64.zip
  unpack: true
```

Then as soon as your app is executed, let `bindeps` check for and install any missing dependencies.

```ruby
require 'bindeps'

Bindeps.require 'binary_dependencies.yaml'
```

`bindeps` will check run the `version:command` for each dependency. If the return value of the command doesn't match a regular expression test against the `version:number` field, `bindeps` will download the file that matches your system architecture, unpack it unless `unpack` is set to false, and place the binary in your path. Specifically, it is added to the `bin` directory of your RubyGems installation. This means the binary will be in the PATH as long as this version of RubyGems is in use (which is ideal for gem dependencies). If the return value does match, `bindeps` will do nothing.

### Specifying an install directory

Simply pass the destination directory as the second argument to `Bindeps::require()`:

```ruby
Bindeps.require('deps.yaml', '/some/install/directory')

```

## Contributing

1. Fork it
2. Create your feature branch (`git checkout -b my-new-feature`)
3. Commit your changes *with tests* (`git commit -am 'Add some feature'`)
4. Push to the branch (`git push origin my-new-feature`)
5. Create new Pull Request
