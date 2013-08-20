require 'test/unit'
require_relative '../lib/Sequence'
require_relative '../lib/Alignment'
require_relative '../lib/NCBIBlastResult'
require_relative '../lib/AccessionNumGroup'

class AccessionNumGroupTest < Test::Unit::TestCase
    def testCreation
        #both params valid
        acc_num = "XM_512.4"
        expr_sig_len = 8
        actual = AccessionNumGroup.new(acc_num,expr_sig_len)
        assert_not_nil(actual)

        #acc_num invalid
        acc_num = "XM_512.4 other strings"
        expr_sig_len = 8
        assert_raise ArgumentError do
            AccessionNumGroup.new(acc_num,expr_sig_len)
        end

        #expr_sig_len invalid
        acc_num = "XM_512.4"
        expr_sig_len = "string"
        assert_raise ArgumentError do
            AccessionNumGroup.new(acc_num,expr_sig_len)
        end

        #both params invalid
        acc_num = "XM_512.4 other strings"
        expr_sig_len = "string"
        assert_raise ArgumentError do
            AccessionNumGroup.new(acc_num,expr_sig_len)
        end
    end
end
