require 'test/unit'
require_relative '../lib/Sequence'

class SequenceTest < Test::Unit::TestCase
    def testCreation
        #both params valid
        id = "abc:123"
        bp_list = "atcgATCG"
        actual = Sequence.new(id,bp_list)
        assert_not_nil(actual)

        id = "abc:123"
        bp_list = "atcgATCG\n"
        actual = Sequence.new(id,bp_list)
        assert_not_nil(actual)

        #id invalid
        id = "abc 123"
        bp_list = "atcgATCG"
        assert_raise ArgumentError do
            Sequence.new(id,bp_list)
        end

        #bp_list invalid
        id = "abc:123"
        bp_list = "atcgA TCG"
        assert_raise ArgumentError do
            Sequence.new(id,bp_list)
        end

        #both params invalid
        id = "abc1\t23"
        pb_list = "atcgATCGz"
        assert_raise ArgumentError do
            Sequence.new(id,bp_list)
        end
    end
    def testAccessors
        #both params valid
        id = "abc:123"
        bp_list = "atcgATCG"
        seq = Sequence.new(id,bp_list)
        actual_id = seq.id
        expected_id = "abc:123"
        assert_equal(expected_id,actual_id)

        actual_bp_list = seq.bp_list
        expected_bp_list = "atcgATCG"
        assert_equal(expected_bp_list,actual_bp_list)
    end
end
