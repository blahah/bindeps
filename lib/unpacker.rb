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
require 'open3'

module Unpacker

  class UnrecognizedArchiveError < StandardError; end

  SUPPORTED_FILEEXTS = %w[tar rar zip gz bz tgz bgz tar]

  def self.unpack(file, tmpdir = "/tmp", &block)
    Dir.mktmpdir 'unpacker' do |dir|
      cmd = case file
            when /rar$/
              "unrar x -y #{file} #{dir}"
            when /(tar|tgz|tar\.gz|tar\.bz|tbz)$/
              "tar xvf #{file} --directory #{dir}"
            when /zip$/
              "unzip #{file} -d #{dir}"
            when /gz$/
              "gunzip -c #{file} #{File.join(dir, "gz-contents")}"
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
    SUPPORTED_FILEEXTS.include? File.extname(file_name).sub('.', '')
  end

  def self.valid?(file_path, file_name = file_path)
    cmd = case file_name
          when /rar$/
            "unrar t #{file_path}"
          when /(tar|tar\.bz|tbz)$/
            "tar tf #{file_path}"
          when /zip$/
            "zip -T #{file_path}"
          when /gz|tgz$/
            "gunzip -t #{file_path}"
          else
            raise UnrecognizedArchiveError
          end
    `#{cmd}`
    true
  rescue
    false
  end

end # module Unpacker
