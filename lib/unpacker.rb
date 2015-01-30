# Lifted from https://github.com/underlog/unpacker with the reliance on ShellShot removed
#
# Copyright (c) 2009 Petyo Ivanov
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

require 'fileutils'
require 'tmpdir'
require 'fixwhich'
require 'open3'
include Which

module Unpacker

  class UnrecognizedArchiveError < StandardError; end

  def self.unpack(file, tmpdir = "/tmp", &block)
    Dir.mktmpdir 'unpacker' do |dir|
      cmd = case file
            when /(tar|tgz|tar\.gz)$/
              "tar xzf #{file} --directory #{dir}"
            when /(tar\.bz|tbz|tar\.bz2)$/
              "tar xjf #{file} --directory #{dir}"
            when /zip$/
              "unzip #{file} -d #{dir}"
            when /gz$/
              "gunzip #{file}"
            when /bz2$/
              "bunzip #{file}"
            else
              raise UnrecognizedArchiveError
            end
      stdout, stderr, status = Open3.capture3 cmd
      if !status.success?
        raise RuntimeError.new("There was a problem unpacking #{file}\n#{stderr}")
      end
      block.call Dir.new(dir)
    end
  end

  def self.archive?(file_name)
    supported = Unpacker.supported_archives
    ext = File.extname(file_name).sub('.', '')
    return true if ext==""
    support = supported.include? ext
    if !support
      help = case file_name
      when /(tar|tgz|tar\.gz|tar\.bz|tbz|tar\.bz2)$/
        "Please install tar"
      when /zip$/
        "Please install unzip"
      when /gz$/
        "Please install gunzip"
      when /bz2$/
        "Please install bunzip2"
      else
        raise UnrecognizedArchiveError
      end
      msg = "Archive type not supported: #{ext}\n#{help}"
      raise UnrecognizedArchiveError.new(msg)
    end
    support
  end

  def self.valid?(file_path, file_name = file_path)
    cmd = case file_name
          when /(tar|tar\.bz|tbz)$/
            "tar tf #{file_path}"
          when /zip$/
            "zip -T #{file_path}"
          when /gz|tgz$/
            "gunzip -t #{file_path}"
          when /bz2$/
            "bunzip2 -t #{file_path}"
          else
            raise UnrecognizedArchiveError
          end
    stdout, stderr, status = Open3.capture3 cmd
    if status.success?
      true
    else
      false
    end
  end

  def self.supported_archives
    supported = []
    if !which('unrar').empty?
      supported << "rar"
    end
    if !which('tar').empty?
      %w[tar tgz tgz tar.gz tar.bz tar.bz2 tbz].each do |ext|
        supported << ext
      end
    end
    if !which('unzip').empty?
      supported << "zip"
    end
    if !which('gunzip').empty?
      supported << "gz"
    end
    if !which('bunzip2').empty?
      supported << "bz2"
    end
    supported
  end

end # module Unpacker
