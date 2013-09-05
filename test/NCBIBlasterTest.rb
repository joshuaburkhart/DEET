require 'test/unit'
require_relative '../lib/Sequence'
require_relative '../lib/Alignment'
require_relative '../lib/NCBIBlastResult'
require_relative '../lib/NCBIBlaster'

class NCBIBlasterTest < Test::Unit::TestCase
    def setup
        @valid_text_result = <<EOS
Query ID=F5BTJ3O01A05UQ
 
Length=507


Score     E
Sequences producing significant alignments:                       (Bits)  Value  N

ref|XM_309606.5|  Anopheles gambiae str. PEST AGAP004052-PA (A...  95.2    3e-36  4
ref|NM_001170892.1|  Nasonia vitripennis homeobox protein pros...  95.6    2e-23  3
ref|XM_966571.2|  PREDICTED: Tribolium castaneum similar to ho...    101   7e-23  3
gb|APGK01052080.1|  Dendroctonus ponderosae Seq01052090, whole...  96.6    6e-22  3
ref|XM_002427623.1|  Pediculus humanus corporis homeobox prote...  92.0    3e-20  2
EOS
        @invalid_text_result = <<EOS
Query ID=F5BTJ3O01A019K
 
Length=90


<b>No significant similarity found.</b> For reasons why, <A HREF = "Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=FAQ"><b>click here</A>.</b><br><br>
 Database: Nucleotide collection (nt)
   Posted date:  Jul 29, 2013  6:14 AM
 Number of letters in database: 48,173,360,552
 Number of sequences in database:  19,024,455

Lambda      K        H
  0.318    0.134    0.401 
Matrix: BLOSUM62
EOS
        @log_filename = "#{Time.now.to_i}.testlog"
        @loghandl = File.open(@log_filename,"w")
    end
    def teardown
        @loghandl.close
        %x(rm #{@log_filename})
    end
    def testBuildBlastResult
        #both params valid
        seq = Sequence.new("valid-id","atcgATCG")
        text_result = @valid_text_result
        blaster = NCBIBlaster.new(@loghandl)
        actual = blaster.buildNCBIResult(text_result,seq)
        assert_not_nil(actual)
        expected = 2
        assert_equal(expected,actual.alignmentCount)
        assert(actual.valid == true)
        assert_equal(seq,actual.sequence)
        assert_equal("XM_309606.5",actual.bestAlignment.accession_num)

        #result invalid
        seq = Sequence.new("valid-id","atcgATCG")
        text_result = @invalid_text_result
        blaster = NCBIBlaster.new(@loghandl)
        actual = blaster.buildNCBIResult(text_result,seq)
        assert(actual.valid == false)

        #sequence invalid
        seq = nil
        text_result = @valid_text_result
        blaster = NCBIBlaster.new(@loghandl)
        assert_raise ArgumentError do
            blaster.buildNCBIResult(text_result,seq)
        end

        #both params invalid
        seq = nil
        text_result = @invalid_text_result
        blaster = NCBIBlaster.new(@loghandl)
        assert_raise ArgumentError do
            blaster.buildNCBIResult(text_result,seq)
        end
    end
    def testSubmitTblastxQuery
        seq = Sequence.new("valid-id","atctagagagagagtttc")
        blaster = NCBIBlaster.new(@loghandl)
        sleep(1) #don't overload ncbi servers
        actual = blaster.submitTblastxQuery(seq)
        #sample return val (rid will vary)
        #<struct NCBIBlaster::PutResponse rid="14WWC21M014", seq=#<Sequence:0x8d5ea64 @id="valid-id", @bp_list="atctagagagagagtttc">>
        assert_equal(NCBIBlaster::PutResponse,actual.class)
        assert_equal(seq.id,actual.seq.id)
        assert_equal(seq.bp_list,actual.seq.bp_list)
        assert(actual.rid.match(/[0-9A-Z-]+/))
    end
    def testFetchTblastxResult
        #put_response valid
        seq = Sequence.new("valid-id","atctagagagagagtttc")
        blaster = NCBIBlaster.new(@loghandl)
        sleep(1) #don't overload ncbi servers
        put_response = blaster.submitTblastxQuery(seq)
        actual = blaster.fetchTblastxResult(put_response)
        assert_not_nil(actual)
        assert_equal(NCBIBlastResult,actual.class)
        assert(actual.valid == false)

        #put_response real
        seq = Sequence.new("XM_001663206.1-sample","ATGCTATGCTATTCACTTCTAATACTAATTGATACCACCATCTTCTCTACTTCCAGCAATTCAGCTCAACCAGCGCCGTCTGTAAATAGCCACGACACACCGCAGAGCACGCGAAGTTCGGCGAATCCCAACAACAACCGTCTGCGAAGACTACGTTCTGCAAGTACTCAGTCGAGTTCGGCAACCAATTTGAACAGCGTCAACAACAACAACGCCAACCCGACCACTGGCTCCGGTGCCAACTCGAACAACAGCCGGCATCCCCCTAATCGCCATCGGACATCAGGAGCCACCGGAACAAGCACTGCAGCAAACACGGGTGGGATAATCTCAAGACTGGTTGACAGTAGTAGCCACGCAGCAGCAGCAGGAGCAGCGCCAGCAACGACGACGACGACTACGACCGGAAGCAACAATAACCGTGGGTCGCCGGGCCCACAGCGTGAATCTACCCCCGTTCGGAAGAAGG")
        blaster = NCBIBlaster.new(@loghandl)
        sleep(1) #don't overload ncbi servers
        put_response = blaster.submitTblastxQuery(seq)
        actual = blaster.fetchTblastxResult(put_response)
        assert_not_nil(actual)
        assert_equal(NCBIBlastResult,actual.class)
        assert(actual.valid == true)

        #put_response invalid
        put_response = nil
        blaster = NCBIBlaster.new(@loghandl)
        assert_raise ArgumentError do
            blaster.fetchTblastxResult(put_response)
        end
    end
end
