require "bindeps/version"
require "unpacker"
require "fixwhich"
require "tmpdir"
require "yaml"

module Bindeps

  class DownloadFailedError < StandardError; end
  class UnsupportedSystemError < StandardError; end

  def self.require dependencies
    if dependencies.is_a? String
      dependencies = YAML.load_file dependencies
    end
    Dir.mktmpdir do |tmpdir|
      Dir.chdir(tmpdir) do
        dependencies.each_pair do |name, config|
          unpack = config.key?('unpack') ? config['unpack'] : true;
          d = Dependency.new(name,
                             config['binaries'],
                             config['version'],
                             config['url'],
                             unpack)
          d.install_missing
        end
      end
    end
  end

  # Check whether all dependencies are installed. Return an array of missing
  # dependencies.
  def self.missing dependencies
    if dependencies.is_a? String
      dependencies = YAML.load_file dependencies
    end
    missing = []
    Dir.mktmpdir do |tmpdir|
      Dir.chdir(tmpdir) do
        dependencies.each_pair do |name, config|
          unpack = config.key?('unpack') ? config['unpack'] : true;
          d = Dependency.new(name,
                             config['binaries'],
                             config['version'],
                             config['url'],
                             unpack)
          missing << d unless d.all_installed?
        end
      end
    end
    missing
  end

  class Dependency

    attr_reader :name, :version, :binaries

    include Which

    def initialize(name, binaries, versionconfig, urlconfig, unpack)
      @name = name
      unless binaries.is_a? Array
        raise ArgumentError,
              "binaries must be an array"
      end
      @binaries = binaries
      @version = versionconfig['number']
      @version_cmd = versionconfig['command']
      @url = choose_url urlconfig
      @unpack = unpack
    end

    def install_missing
      unless all_installed?
        puts "Installing #{@name} (#{@version})..."
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
      wget = which('wget').first
      curl = which('curl').first
      if wget
        cmd = "#{wget} #{@url}"
        stdout, stderr, status = Open3.capture3 cmd
      elsif curl
        cmd = "#{curl} -O -J -L #{@url}"
        stdout, stderr, status = Open3.capture3 cmd
      else
        msg = "You don't have curl or wget?! What kind of computer is "
        msg << "this?! Windows?! BeOS? OS/2?"
        raise DownloadFailedError.new(msg)
      end
      if !status.success?
        raise DownloadFailedError,
              "download of #{@url} for #{@name} failed:\n#{stdout}\n#{stderr}"
      end
    end

    def unpack
      archive = File.basename(@url)
      Unpacker.archive? archive
      if @unpack
        Unpacker.unpack(archive) do |dir|
          Dir.chdir dir do
            Dir['**/*'].each do |extracted|
              if @binaries.include? File.basename(extracted)
                install(extracted) unless File.directory?(extracted)
              end
            end
          end
        end
      else
        install(@binaries.first)
      end
    end

    def all_installed?
      @binaries.each do |bin|
        unless installed? bin
          return false
        end
      end
      true
    end

    def installed? bin
      path = which(bin)
      if path.length > 0
        ret = `#{@version_cmd} 2>&1`.split("\n").map{ |l| l.strip }.join('|')
        if ret && (/#{@version}/ =~ ret)
          return path
        end
      else
        if Dir.exist?("#{ENV['HOME']}/.local/bin")
          ENV['PATH'] = ENV['PATH'] + ":#{ENV['HOME']}/.local/bin"
          path = which(bin)
          if path.length > 0
            return path
          end
        end
      end
      false
    end

    def install bin
      gem_home = ENV['GEM_HOME']
      home = ENV['HOME']
      bindir = "#{home}/.local/bin"
      if gem_home.nil?
        if !Dir.exist?("#{home}/.local")
          Dir.mkdir("#{home}/.local")
        end
        if !Dir.exist?("#{home}/.local/bin")
          Dir.mkdir("#{home}/.local/bin")
        end
        ENV['PATH'] = ENV['PATH'] + ":#{ENV['HOME']}/.local/bin"
      else
        bindir = File.join(ENV['GEM_HOME'], 'bin')
      end
      install_location = File.join(bindir, File.basename(bin))
      FileUtils.install(bin, install_location, :mode => 0775)
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
