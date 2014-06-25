require "bindeps/version"
require "unpacker"
require "which"
require "tmpdir"
require "yaml"

module Bindeps

  class UnsupportedArchiveError < StandardError; end
  class DownloadFailedError < StandardError; end
  class UnsupportedSystemError < StandardError; end

  def self.require dependencies
    if dependencies.is_a? String
      dependencies = YAML.load_file dependencies
    end
    Dir.mktmpdir do |tmpdir|
      Dir.chdir(tmpdir) do
        dependencies.each_pair do |name, config|
          d = Dependency.new(name,
                             config[:binaries],
                             config[:version],
                             config[:url])
          d.install_missing
        end
      end
    end
  end

  class Dependency

    def init(name, binaries, versionconfig, urlconfig)
      @name = name
      @binaries = binaries
      @version = versionconfig[:number]
      @version_cmd = versionconfig[:command]
      @url = resolve_url urlconfig
    end

    def install_missing
      unless all_installed?
        binaries.each do |bin|
            download
            unpack
            return
          end
        end
      end
    end

    def choose_url urlconfig
      if urlconfig.key? System.arch
        sys = urlconfig[System.arch]
        if sys.key? System.arch
          return sys[System.arch]
        else
          raise UnsupportedSystemError,
                "bindeps config for #{@name} doesn't contain an entry for #{System.arch}"
        end
      else
        raise UnsupportedSystemError,
              "bindeps config for #{@name} doesn't contain an entry for #{System.os}"
      end
    end

    def download
      `curl -O -J -L #{@url}`
      unless $?.to_i == 0
        raise DownloadFailedError,
              "download of #{@url} for #{@name} failed"
    end

    def unpack
      Unpacker.unpack(archive) do |dir|
        Dir.chdir dir do
          Dir['*'].each do |extracted|
            if @binaries.include? extracted
              install extracted
            end
          end
        end
      end
    end

    def all_installed?
      binaries.each do |bin|
        return false unless installed? bin
      end
      true
    end

    def installed? bin
      `#{version_cmd}` Which.which bin
    end

    def install bin
      bindir = File.join(ENV['GEM_HOME'], 'bin')
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
