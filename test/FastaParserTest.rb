require 'test/unit'
require_relative '../lib/FastaParser'
require_relative '../lib/Sequence'

class SequenceTest < Test::Unit::TestCase
    V_FASTA_N = "valid_fasta.fasta"
    I_FASTA_N = "invalid_fasta.fasta"
    def setup
        valid_fasta_contents = <<EOS
>sample_seq_id
ATCTATCG
TTTGGGA
TTTTTTTT
>sample_seq_id2
atcATCaaatttatatatatatatatat
EOS
        valid_fasta_handl = File.open(V_FASTA_N,"w")
        valid_fasta_handl.puts(valid_fasta_contents)
        valid_fasta_handl.close

        invalid_fasta_contents = <<EOS
>sam ple seq_id
AZFG
>
\n
>valid_id
ZZZZMMMM
>in valid id
ATCGATCG
>VALID
ATCGATCGATCT
EOS
        invalid_fasta_handl = File.open(I_FASTA_N,"w")
        invalid_fasta_handl.puts(invalid_fasta_contents)
        invalid_fasta_handl.close
        @log_filename = "#{Time.now.to_i}.testlog"
        @loghandl = File.open(@log_filename,"w")
    end
    def teardown
        %x(rm -f #{V_FASTA_N})
        %x(rm -f #{I_FASTA_N})
        @loghandl.close
        %x(rm #{@log_filename})
    end
    def testCreation
        #both params valid
        filename = V_FASTA_N
        min_len = 20   
        actual = FastaParser.new(filename,min_len,@loghandl)
        assert_not_nil(actual)
        assert_equal(filename,actual.fasta_filename)
        assert_equal(min_len,actual.min_len)

        filename = I_FASTA_N
        min_len = 20
        actual = FastaParser.new(filename,min_len,@loghandl)
        assert_not_nil(actual)
        assert_equal(filename,actual.fasta_filename)
        assert_equal(min_len,actual.min_len)

        #filename invalid
        filename = "nonexistant file"
        min_len = 20
        assert_raise ArgumentError do
            FastaParser.new(filename,min_len,@loghandl)
        end

        #min_len invalid
        filename = I_FASTA_N
        min_len = 'a'
        assert_raise ArgumentError do
            FastaParser.new(filename,min_len,@loghandl)
        end

        filename = V_FASTA_N
        min_len = "2.91232354"
        assert_raise ArgumentError do
            FastaParser.new(filename,min_len,@loghandl)
        end

        #both params invalid
        filename = "another-nonexistant, , file"
        min_len = "3.141592653589"
        assert_raise ArgumentError do
            FastaParser.new(filename,min_len,@loghandl)
        end
    end
    def testFileOpenClose
        filename = V_FASTA_N
        min_len = 20   
        parser = FastaParser.new(filename,min_len,@loghandl)
        parser.open
        actual = parser.fasta_filehandl
        expected = File.open(V_FASTA_N,"w")
        assert_not_nil(parser.fasta_filehandl)
        assert_equal(expected.inspect,actual.inspect)
        parser.close
        expected.close
        assert_equal(expected.inspect,actual.inspect)
    end
    def testNextSeq
        #valid fasta
        filename = V_FASTA_N
        min_len = 8
        parser = FastaParser.new(filename,min_len,@loghandl)
        parser.open

        actual_seq1 = parser.nextSeq
        expected_seq1 = Sequence.new("sample_seq_id","ATCTATCGTTTGGGATTTTTTTT")
        assert_equal(Sequence,actual_seq1.class)
        assert_equal(expected_seq1.id,actual_seq1.id)
        assert_equal(expected_seq1.bp_list,actual_seq1.bp_list)

        actual_seq2 = parser.nextSeq
        expected_seq2 = Sequence.new("sample_seq_id2","atcATCaaatttatatatatatatatat")
        assert_equal(Sequence,actual_seq2.class)
        assert_equal(expected_seq2.id,actual_seq2.id)
        assert_equal(expected_seq2.bp_list,actual_seq2.bp_list)

        actual_seq3 = parser.nextSeq
        expected_seq3 = nil
        assert_equal(expected_seq3,actual_seq3)
        parser.close

        actual_seq4 = parser.nextSeq
        expected_seq4 = nil
        assert_equal(expected_seq4,actual_seq4)

        #invalid fasta
        filename = I_FASTA_N
        min_len = 8
        parser = FastaParser.new(filename,min_len,@loghandl)
        parser.open

        assert_raise ArgumentError do
            actual_seq1 = parser.nextSeq
        end

        parser.close

        #increasing min_len
        filename = V_FASTA_N
        min_len = 9
        parser = FastaParser.new(filename,min_len,@loghandl)
        parser.open

        parser.nextSeq

        actual_seq1 = parser.nextSeq
        expected_seq1 = Sequence.new("sample_seq_id2","atcATCaaatttatatatatatatatat")
        assert_equal(Sequence,actual_seq1.class)
        assert_equal(expected_seq1.id,actual_seq1.id)
        assert_equal(expected_seq1.bp_list,actual_seq1.bp_list)

        actual_seq2 = parser.nextSeq
        expected_seq2 = nil
        assert_equal(expected_seq2,actual_seq2)
        parser.close

        actual_seq3 = parser.nextSeq
        expected_seq3 = nil
        assert_equal(expected_seq3,actual_seq3)
    end
end
