class SimpleFastaParser
    @fasta_filename
    @fasta_filehandl
    @min_len
    def initialize(fasta_filename,min_len)
        @fasta_filename = fasta_filename
        @min_len = min_len
    end
    def open()
        @fasta_filehandl = File.open(fasta_filename,"r")
    end
    def nextSeq()
        seq_id = nil
        seq_bp_list = ""
        seq = nil
        begin        
            line1 = @fasta_filehandl.gets
            if(line1.match(/^>(\w*)/)
               seq_id = $1
               line2 = @fasta_filehandl.gets
               if(line2.match(/^([atcgATCG]+)/)
                  seq_bp_list = $1
               end
            end
        end while(seq_bp_list.length < @min_len && !@fasta_filehandl.eof)
        seq = Sequence.new(seq_id,seq_bp_list)
        return seq
    end
    def close()
        @fasta_filehandl.close
    end
end
