require "bindeps/version"
require "unpacker"
require "which_works"
require "tmpdir"
require "yaml"

module Bindeps

  class DownloadFailedError < StandardError; end
  class UnsupportedSystemError < StandardError; end

  def self.require dependencies
    if dependencies.is_a? String
      dependencies = YAML.load_file dependencies
    end
    tmpdir = Dir.mktmpdir
    Dir.chdir(tmpdir) do
      dependencies.each_pair do |name, config|
        d = Dependency.new(name,
                           config['binaries'],
                           config['version'],
                           config['url'])
        d.install_missing
      end
    end
  end

  class Dependency

    def initialize(name, binaries, versionconfig, urlconfig)
      @name = name
      unless binaries.is_a? Array
        raise ArgumentError,
              "binaries must be an array"
      end
      @binaries = binaries
      @version = versionconfig['number']
      @version_cmd = versionconfig['command']
      @url = choose_url urlconfig
    end

    def install_missing
      unless all_installed?
        puts "binary dependency #{@name} not installed"
        puts "it will now be downloaded and installed"
        download
        unpack
      end
    end

    def choose_url urlconfig
      arch = System.arch
      if urlconfig.key? arch
        sys = urlconfig[arch]
        os = System.os.to_s
        if sys.key? os
          return sys[os]
        else
          raise UnsupportedSystemError,
                "bindeps config for #{@name} #{arch} doesn't contain an entry for #{os}"
        end
      else
        raise UnsupportedSystemError,
              "bindeps config for #{@name} doesn't contain an entry for #{arch}"
      end
    end

    def download
      `curl -O -J -L #{@url}`
      unless $?.to_i == 0
        raise DownloadFailedError,
              "download of #{@url} for #{@name} failed"
      end
    end

    def unpack
      archive = File.basename(@url)
      Unpacker.unpack(archive) do |dir|
        Dir.chdir dir do
          Dir['**/*'].each do |extracted|
            if @binaries.include? File.basename(extracted)
              install(extracted) unless File.diretory?(extracted)
            end
          end
        end
      end
    end

    def all_installed?
      @binaries.each do |bin|
        unless installed? bin
          puts "required binary #{bin} is not installed"
          return false
        end
      end
      true
    end

    def installed? bin
      path = Which.which(bin)
      if path
        ret = `#{@version_cmd}`.strip
        if ret && (/#{@version}/ =~ ret)
          return path
        end
      end
      false
    end

    def install bin
      bindir = File.join(ENV['GEM_HOME'], 'bin')
      puts "installing #{bin} to #{bindir}"
      FileUtils.cp(bin, File.join(bindir, File.basename(bin)))
    end

  end

  class System

    require 'rbconfig'

    def self.os
      (
        host_os = RbConfig::CONFIG['host_os']
        case host_os
        when /mswin|msys|mingw|cygwin|bccwin|wince|emc/
          :windows
        when /darwin|mac os/
          :macosx
        when /linux/
          :linux
        when /solaris|bsd/
          :unix
        else
          raise UnsupportedSystemError,
                "can't install #{@name}, unknown os: #{host_os.inspect}"
        end
      )
    end

    def self.arch
      Gem::Platform.local.cpu == 'x86_64' ? '64bit' : '32bit'
    end

  end

end
