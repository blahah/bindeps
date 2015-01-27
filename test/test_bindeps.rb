require 'helper'

class TestBindeps < Test::Unit::TestCase

  context "bindeps" do

    setup do
      test_dir = File.dirname(__FILE__)
      @data_dir = File.join(test_dir, 'data')
    end

    teardown do
      # delete fake binaries from
      bindir = File.join(ENV['GEM_HOME'], 'bin')
      `rm #{bindir}/fakebin` if File.exist?("#{bindir}/fakebin")
      `rm #{bindir}/fakebin2` if File.exist?("#{bindir}/fakebin2")
      `rm #{bindir}/fakebin3` if File.exist?("#{bindir}/fakebin3")
      `rm #{bindir}/fakebin4` if File.exist?("#{bindir}/fakebin4")
    end

    should "identify and install missing dependencies" do
      test_yaml = File.join(@data_dir, 'fakebin.yaml')
      Bindeps.require test_yaml
      assert_equal 'success', `fakebin`.strip
    end

    should "identify dependencies that are missing" do
      test_yaml = File.join(@data_dir, 'neverinstalled.yaml')
      missing = Bindeps.missing test_yaml
      assert_equal 1, missing.length
      assert_equal 'neverinstalled', missing.first.binaries.first
    end

    should "handle case where version is not on first line" do
      test_yaml = File.join(@data_dir, 'fakebin2.yaml')
      # install fakebin2
      Bindeps.require test_yaml
      # now Dependency should detect it as installed
      deps = YAML.load_file test_yaml
      deps.each_pair do |name, config|
        d = Bindeps::Dependency.new(name,
                                    config['binaries'],
                                    config['version'],
                                    config['url'],
                                    config['unpack'])
        assert d.installed?('fakebin2')
      end
    end

    should "handle version output to stderr" do
      test_yaml = File.join(@data_dir, 'fakebin3.yaml')
      # install fakebin3
      Bindeps.require test_yaml
      # now Dependency should detect it as installed
      deps = YAML.load_file test_yaml
      deps.each_pair do |name, config|
        d = Bindeps::Dependency.new(name,
                                    config['binaries'],
                                    config['version'],
                                    config['url'],
                                    config['unpack'])
        assert d.installed?('fakebin3')
        assert d.all_installed?
      end
    end

    should "fail when no download is specified for the local system" do
      require 'rbconfig'
      original_os = RbConfig::CONFIG['host_os']
      RbConfig::CONFIG['host_os'] = 'fail'
      assert_raise Bindeps::UnsupportedSystemError do
        Bindeps::System.os
      end
      urlconfig = { '64bit' => { 'unix' => 'url' } }
      versionconfig = { 'number' => 1, 'command' => 'getversion' }
      RbConfig::CONFIG['host_os'] = 'linux'
      assert_raise Bindeps::UnsupportedSystemError do
        Bindeps::Dependency.new('test', ['binary'],
                                versionconfig, urlconfig, false)
      end
      RbConfig::CONFIG['host_os'] = 'solaris'
      original_cpu = Gem::Platform.local.cpu
      Gem::Platform.local.cpu = 'fail'
      assert_raise Bindeps::UnsupportedSystemError do
        Bindeps::Dependency.new('test', ['binary'],
                                versionconfig, urlconfig, false)
      end
      # set os info back to the truth
      RbConfig::CONFIG['host_os'] = original_os
      Gem::Platform.local.cpu = original_cpu
    end

    should "initialize" do
      assert_nothing_raised do
        Bindeps::Dependency.new(
          'test',                     # name
          ['binary'],                 # binaries
          {                           # versionconfig
            'number' => 1,
            'command' => 'getversion'
          },
          {                           # urlconfig
            '64bit' => {
              'linux' => 'url',
              'unix' => 'url',
              'macosx' => 'url',
              'windows' => 'url',
            }
          },
          false                       # unpack
        )
      end
      assert_raise ArgumentError do
        Bindeps::Dependency.new(
          'test',                     # name
          'binary',                   # binaries is no longer an array
          {                           # versionconfig
            'number' => 1,
            'command' => 'getversion'
          },
          {                           # urlconfig
            '64bit' => {
              'linux' => 'url',
              'unix' => 'url',
              'macosx' => 'url',
              'windows' => 'url',
            }
          },
          false                       # unpack
        )
      end
    end

    should "handle binaries that don't need to be unpacked" do
      test_yaml = File.join(@data_dir, 'fakebin4.yaml')
      Bindeps.require test_yaml
      deps = YAML.load_file test_yaml
      deps.each_pair do |name, config|
        d = Bindeps::Dependency.new(name,
                                    config['binaries'],
                                    config['version'],
                                    config['url'],
                                    config['unpack'])
        assert d.installed?('fakebin4'), "fakebin4 installed"
      end
    end

  end

end
