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
    def testAddRes
        #single sequence, single signature, no alignments
        seq = Sequence.new("valid-id","ATCGATCG")
        res = NCBIBlastResult.new(seq)
        expr_sig = "100111010"
        acc_num = "AG_123.4"
        acc_num_grp = AccessionNumGroup.new(acc_num,expr_sig.length)
        assert_raise ArgumentError do
            acc_num_grp.addRes(expr_sig,res)
        end

        #single sequence, single signature, two alignments / sequence
        seq = Sequence.new("valid-id","ATCGATCG")
        res = NCBIBlastResult.new(seq)

        accession_num1 = "DQ083361.1"
        e_value1 = "2.5e-10"
        alignment1 = Alignment.new(seq,accession_num1,e_value1)
        res.addAlignment(alignment1)

        accession_num2 = "X5234234.9"
        e_value2 = "1"
        alignment2 = Alignment.new(seq,accession_num2, e_value2)
        res.addAlignment(alignment2)

        expr_sig = "100111010"
        acc_num_grp = AccessionNumGroup.new(res.bestAlignment.accession_num,expr_sig.length)
        acc_num_grp.addRes(expr_sig,res)
        assert_equal(seq,acc_num_grp.getRepresentativeSeq(expr_sig))
        assert_equal(0,acc_num_grp.getParalogExprSigs.length)
        assert_equal(expr_sig,acc_num_grp.getGeneExprSig)

        #two sequences, single signature, two alignments / sequence
        seq1 = Sequence.new("valid-id_1","ATCGATCGGGGGGGGGG")
        res1 = NCBIBlastResult.new(seq1)

        accession_num1 = "B5204200.4"
        e_value1 = "2.5e-10"
        alignment1 = Alignment.new(seq1,accession_num1,e_value1)
        res1.addAlignment(alignment1)

        accession_num2 = "X5234234.9"
        e_value2 = "1"
        alignment2 = Alignment.new(seq1,accession_num2, e_value2)
        res1.addAlignment(alignment2)

        seq2 = Sequence.new("valid-id_2","ATCGATCGGTGT")
        res2 = NCBIBlastResult.new(seq2)

        accession_num3 = "Z1083361.3"
        e_value3 = "2.5e-8"
        alignment3 = Alignment.new(seq2,accession_num3,e_value3)
        res2.addAlignment(alignment3)

        accession_num4 = "B5204200.4"
        e_value4 = "1.5e-19"
        alignment4 = Alignment.new(seq2,accession_num4, e_value4)
        res2.addAlignment(alignment4)

        expr_sig = "100111010"
        acc_num_grp = AccessionNumGroup.new(res1.bestAlignment.accession_num,expr_sig.length)
        acc_num_grp.addRes(expr_sig,res1)
        acc_num_grp.addRes(expr_sig,res2)
        assert_equal(seq1,acc_num_grp.getRepresentativeSeq(expr_sig))
        assert_equal(0,acc_num_grp.getParalogExprSigs.length)
        assert_equal(expr_sig,acc_num_grp.getGeneExprSig)

        #four sequences, two signatures, two alignments / sequence
        seq1 = Sequence.new("valid-id_1","ATCGATCGGGGGGGGGG")
        res1 = NCBIBlastResult.new(seq1)
        accession_num1 = "B5204200.4"
        e_value1 = "2.5e-10"
        alignment1 = Alignment.new(seq1,accession_num1,e_value1)
        res1.addAlignment(alignment1)
        accession_num2 = "X5234234.9"
        e_value2 = "1"
        alignment2 = Alignment.new(seq1,accession_num2, e_value2)
        res1.addAlignment(alignment2)

        seq2 = Sequence.new("valid-id_2","ATCGATCGGTGT")
        res2 = NCBIBlastResult.new(seq2)
        accession_num3 = "Z1083361.3"
        e_value3 = "2.5e-8"
        alignment3 = Alignment.new(seq2,accession_num3,e_value3)
        res2.addAlignment(alignment3)
        accession_num4 = "B5204200.4"
        e_value4 = "1.5e-19"
        alignment4 = Alignment.new(seq2,accession_num4, e_value4)
        res2.addAlignment(alignment4)

        expr_sig1 = "100111010"

        seq3 = Sequence.new("valid-id_3","ATCGATCGGGGGGATCGGGGG")
        res3 = NCBIBlastResult.new(seq3)
        accession_num5 = "B5204200.4"
        e_value5 = "2.5e-100"
        alignment5 = Alignment.new(seq3,accession_num5,e_value5)
        res3.addAlignment(alignment5)
        accession_num6 = "XX234234.1"
        e_value6 = "1.3e-2"
        alignment6 = Alignment.new(seq3,accession_num6, e_value6)
        res3.addAlignment(alignment6)

        seq4 = Sequence.new("valid-id_4","ATTAAATCGATCGGTGT")
        res4 = NCBIBlastResult.new(seq4)
        accession_num7 = "B1083361.3"
        e_value7 = "2.5e-18"
        alignment7 = Alignment.new(seq2,accession_num7,e_value7)
        res4.addAlignment(alignment7)
        accession_num8 = "B5204200.4"
        e_value8 = "1.5e-190"
        alignment8 = Alignment.new(seq2,accession_num8, e_value8)
        res4.addAlignment(alignment8)

        expr_sig2 = "111111010"
        acc_num_grp = AccessionNumGroup.new(res1.bestAlignment.accession_num,expr_sig1.length)
        acc_num_grp.addRes(expr_sig1,res1)
        acc_num_grp.addRes(expr_sig1,res2)
        acc_num_grp.addRes(expr_sig2,res3)
        acc_num_grp.addRes(expr_sig2,res4)
        assert_equal(seq1,acc_num_grp.getRepresentativeSeq(expr_sig1))
        assert_equal(seq3,acc_num_grp.getRepresentativeSeq(expr_sig2))
        assert_equal(1,acc_num_grp.getParalogExprSigs.length)
        assert_equal(expr_sig2,acc_num_grp.getGeneExprSig)

        #four sequences, two signatures, 1.5 alignments / sequence (one sequence with no alignments)
        seq1 = Sequence.new("valid-id_1","ATCGATCGGGGGGGGGG")
        res1 = NCBIBlastResult.new(seq1)
        accession_num1 = "B5204200.4"
        e_value1 = "2.5e-10"
        alignment1 = Alignment.new(seq1,accession_num1,e_value1)
        res1.addAlignment(alignment1)
        accession_num2 = "X5234234.9"
        e_value2 = "1"
        alignment2 = Alignment.new(seq1,accession_num2, e_value2)
        res1.addAlignment(alignment2)

        seq2 = Sequence.new("valid-id_2","ATCGATCGGTGT")
        res2 = NCBIBlastResult.new(seq2)
        accession_num3 = "Z1083361.3"
        e_value3 = "2.5e-8"
        alignment3 = Alignment.new(seq2,accession_num3,e_value3)
        res2.addAlignment(alignment3)
        accession_num4 = "B5204200.4"
        e_value4 = "1.5e-19"
        alignment4 = Alignment.new(seq2,accession_num4, e_value4)
        res2.addAlignment(alignment4)

        expr_sig1 = "100111010"

        seq3 = Sequence.new("valid-id_3","ATCGATCGGGGGGATCGGGGG")
        res3 = NCBIBlastResult.new(seq3)

        seq4 = Sequence.new("valid-id_4","ATTAAATCGATCGGTGT")
        res4 = NCBIBlastResult.new(seq4)
        accession_num7 = "B1083361.3"
        e_value7 = "2.5e-18"
        alignment7 = Alignment.new(seq2,accession_num7,e_value7)
        res4.addAlignment(alignment7)
        accession_num8 = "B5204200.4"
        e_value8 = "1.5e-190"
        alignment8 = Alignment.new(seq2,accession_num8, e_value8)
        res4.addAlignment(alignment8)

        expr_sig2 = "111111010"
        acc_num_grp = AccessionNumGroup.new(res1.bestAlignment.accession_num,expr_sig1.length)
        acc_num_grp.addRes(expr_sig1,res1)
        acc_num_grp.addRes(expr_sig1,res2)
        assert_raise ArgumentError do
            acc_num_grp.addRes(expr_sig2,res3)
        end
        acc_num_grp.addRes(expr_sig2,res4)
        assert_equal(seq1,acc_num_grp.getRepresentativeSeq(expr_sig1))
        assert_equal(seq4,acc_num_grp.getRepresentativeSeq(expr_sig2))
        assert_equal(1,acc_num_grp.getParalogExprSigs.length)
        assert_equal(expr_sig2,acc_num_grp.getGeneExprSig)
    end
end
