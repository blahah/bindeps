require 'helper'

class TestBinDeps < Test::Unit::TestCase

  context "bindeps" do

    setup do
      test_dir = File.dirname(__FILE__)
      @data_dir = File.join(test_dir, 'data')
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
      assert_equal 'neverinstalled', missing.first
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
