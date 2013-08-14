require 'test/unit'
require_relative '../lib/Sequence'
require_relative '../lib/Alignment'
require_relative '../lib/NCBIBlastResult'

class NCBIBlastResultTest < Test::Unit::TestCase
    def testCreation
        #sequence valid
        seq = Sequence.new("valid-id","atcgATCG")
        actual = NCBIBlastResult.new(seq)
        assert_not_nil(actual)

        #sequence invalid
        seq = nil
        assert_raise ArgumentError do
            NCBIBlastResult.new(seq)
        end
    end
    def testAddAlignment
        #single alignment
        seq = Sequence.new("valid-id","atcgATCG")
        align = Alignment.new(seq,"XM_1234.5A","2.3e-5")
        result = NCBIBlastResult.new(seq)
        result.addAlignment(align)
        actual = result.alignments
        expected = [align]
        assert_equal(expected.class,actual.class)
        assert_equal(expected.length,actual.length)
        assert(seqAryCmp(expected,actual))

        #multiple alignments
        seq1 = Sequence.new("valid-id","atcgATCG")
        align1 = Alignment.new(seq,"XM_1234.5A","2.3e-5")
        align2 = Alignment.new(seq,"XM_1234.5B","2.3")
        result = NCBIBlastResult.new(seq)
        result.addAlignment(align1)
        result.addAlignment(align2)
        actual = result.alignments
        expected = [align1,align2]
        assert_equal(expected.class,actual.class)
        assert_equal(expected.length,actual.length)
        assert(seqAryCmp(expected,actual))
    end
    def testAlignCount
        #no alignments
        seq = Sequence.new("valid-id","atcgATCG")
        align = Alignment.new(seq,"XM_1234.5A","2.3e-5")
        result = NCBIBlastResult.new(seq)
        actual = result.alignmentCount
        expected = 0
        assert_equal(expected,actual)

        #single alignment
        seq = Sequence.new("valid-id","atcgATCG")
        align = Alignment.new(seq,"XM_1234.5A","2.3e-5")
        result = NCBIBlastResult.new(seq)
        result.addAlignment(align)
        actual = result.alignmentCount
        expected = 1
        assert_equal(expected,actual)

        #multiple alignments
        seq1 = Sequence.new("valid-id","atcgATCG")
        align1 = Alignment.new(seq,"XM_1234.5A","2.3e-5")
        align2 = Alignment.new(seq,"XM_1234.5B","2.3")
        result = NCBIBlastResult.new(seq)
        result.addAlignment(align1)
        result.addAlignment(align2)
        actual = result.alignmentCount
        expected = 2
        assert_equal(expected,actual)
    end
    def testBestAlignment
        #no alignments
        seq = Sequence.new("valid-id","atcgATCG")
        result = NCBIBlastResult.new(seq)
        actual = result.bestAlignment
        expected = nil
        assert_equal(expected,actual)

        #single alignment
        seq = Sequence.new("valid-id","atcgATCG")
        align = Alignment.new(seq,"XM_1234.5A","2.3e-5")
        result = NCBIBlastResult.new(seq)
        result.addAlignment(align)
        actual = result.bestAlignment
        expected = align
        assert_equal(expected.class,actual.class)
        assert_equal(expected.instance_variables,actual.instance_variables)

        #multiple alignments
        seq1 = Sequence.new("valid-id","atcgATCG")
        align1 = Alignment.new(seq,"XM_1234.5A","2.3e-5")
        align2 = Alignment.new(seq,"XM_1234.5B","2.3")
        result = NCBIBlastResult.new(seq)
        result.addAlignment(align2)
        result.addAlignment(align1)
        actual = result.bestAlignment
        expected = align1
        assert_equal(expected.class,actual.class)
        assert_equal(expected.instance_variables,actual.instance_variables)
    end
    def seqAryCmp(ary1,ary2)
        if(ary1.class == Array && ary2.class == Array)
            ary1.sort! {|a,b| a.e_value <=> b.e_value}
            ary2.sort! {|a,b| a.e_value <=> b.e_value}
            if(ary1.length == ary2.length)
                for i in 0..(ary1.length - 1)
                    if(ary1[i].class != ary2[i].class)
                        return false
                    end
                    if(ary1[i].instance_variables != ary2[i].instance_variables)
                        return false
                    end
                    for j in 0..(ary1[i].instance_variables.length - 1)
                        var1 = ary1[i].instance_variables[j]
                        var2 = ary2[i].instance_variables[j]
                        value1 = ary1[i].instance_variable_get("#{var1}")
                        value2 = ary2[i].instance_variable_get("#{var2}")
                        if(value1 != value2)
                            return false
                        end
                    end
                end
                return true
            end
        end 
        return false
    end
end
