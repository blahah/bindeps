# Bindeps

Simple binary dependency management for Ruby gems

## Installation

`gem install bindeps`

### Using Bundler

Add this line to your application's Gemfile:

    gem 'bindeps'

And then execute:

    $ bundle

Or install it yourself as:

    $ gem install bindeps

## Usage

Create a YAML file describing your dependencies as a dictionary. Read the [bindeps YAML format specifications](wiki/bindeps_YAML_format_specifications).

```yaml
blastplus:
  binaries:
    - makeblastdb
    - blastn
    - tblastn
    - blastp
    - blastx
  version:
    number: '2.2.29'
    command: 'blastx -version'
  url:
    64bit:
      osx: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.29/ncbi-blast-2.2.29+-universal-macosx.tar.gz
      linux: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.29/ncbi-blast-2.2.29+-x64-linux.tar.gz
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
      osx: http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.2.3/bowtie2-2.2.3-macos-x86_64.zip
```

Then as soon as your app is executed, let `bindeps` check for and install any missing dependencies.

```ruby
require 'bindeps'

Bindeps.require 'binary_dependencies.yaml'
```

`bindeps` will check run the `versioncmd` for each dependency. If the return value of the command doesn't match a regular expression test against the `version` field, `bindeps` will download the file that matches your system architecture, unpack it and place the binary in your path. If the return value does match, `bindeps` will do nothing.

## Contributing

1. Fork it
2. Create your feature branch (`git checkout -b my-new-feature`)
3. Commit your changes *with tests* (`git commit -am 'Add some feature'`)
4. Push to the branch (`git push origin my-new-feature`)
5. Create new Pull Request
