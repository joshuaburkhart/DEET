require 'test/unit'
require_relative '../lib/Sequence'
require_relative '../lib/Alignment'

class AlignmentTest < Test::Unit::TestCase
    def testCreation
        #all params valid
        seq = Sequence.new("valid_id","atcgATCG")
        accession_num = "DQ083361"
        e_value = "2.3e-20"
        actual = Alignment.new(seq,accession_num,e_value)
        assert_not_nil(actual)

        seq = Sequence.new("valid_id2","ATATATATatcgATCG")
        accession_num = "DQ_08336.1"
        e_value = "2.3"
        actual = Alignment.new(seq,accession_num,e_value)
        assert_not_nil(actual)

        #sequence invalid
        seq = nil
        accession_num = "DQ_08336.1"
        e_value = "2.3"
        assert_raise ArgumentError do
            Alignment.new(seq,accession_num,e_value)
        end

        seq = "hello1234"
        accession_num = "DQ_08336.1"
        e_value = "2.3"
        assert_raise ArgumentError do
            Alignment.new(seq,accession_num,e_value)
        end

        #accession_num invalid
        seq = Sequence.new("valid_id","atcgATCG")
        accession_num = nil
        e_value = "2.3e-20"
        assert_raise ArgumentError do
            Alignment.new(seq,accession_num,e_value)
        end

        seq = Sequence.new("valid_id","atcgATCG")
        accession_num = 5
        e_value = "2.3e-20"
        assert_raise ArgumentError do
            Alignment.new(seq,accession_num,e_value)
        end

        seq = Sequence.new("valid_id","atcgATCG")
        accession_num = "!@#\$%^&*()"
        e_value = "2.3e-20"
        assert_raise ArgumentError do
            Alignment.new(seq,accession_num,e_value)
        end

        #e_value invalid
        seq = Sequence.new("valid_id","atcgATCG")
        accession_num = "XM_AB893.34A"
        e_value = nil
        assert_raise ArgumentError do
            Alignment.new(seq,accession_num,e_value)
        end

        seq = Sequence.new("valid_id","atcgATCG")
        accession_num = "XM_AB893.34A"
        e_value = 3.45
        assert_raise ArgumentError do
            Alignment.new(seq,accession_num,e_value)
        end

        seq = Sequence.new("valid_id","atcgATCG")
        accession_num = "XM_AB893.34A"
        e_value = "three"
        assert_raise ArgumentError do
            Alignment.new(seq,accession_num,e_value)
        end
    end
end
