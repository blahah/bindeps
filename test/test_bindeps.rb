require 'helper'

class TestBinDeps < Test::Unit::TestCase

  context "bindeps" do

    setup do
      test_dir = File.dirname(__FILE__)
      data_dir = File.join(test_dir, 'data')
      @test_yaml = File.join(data_dir, 'fakebin.yaml')
    end

    # teardown do
    #
    # end
    #
    # should "check if dependencies are installed" do
    #
    # end
    #
    # should "download and unpack dependencies" do
    # end

    should "identify and install missing dependencies" do
      Bindeps.require @test_yaml
      assert_equal 'success', `fakebin`
    end

  end

end
