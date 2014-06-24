require "bindeps/version"
require "unpacker"
require "which"

module Bindeps

  class UnsupportedArchiveError < StandardError; end
  class DownloadFailedError < StandardError; end

  def self.require dependencies
    dependencies.each do |dep|
      Dependency.new(dep[:binaries], dep[:urlconfig])
    end
  end

  class Dependency

    def init(binaries, urlconfig)
      @binaries = binaries
      @url = urlconfig[:url]
    end

    def install_missing
      binaries.each do |bin|
        if !installed? bin
          download
          unpack
          return
        end
      end
    end

    def download
      `curl -O -J -L #{url}`
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
      Which.which bin
    end

    def install bin
      bindir = ENV['GEM_HOME']
      FileUtils.cp(bin, File.join(bindir, File.basename(bin)))
    end

  end


end
