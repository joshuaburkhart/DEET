class SimpleFastaParser
    @fasta_filename
    @fasta_filehandl
    @min_len
    def initialize(fasta_filename,min_len)
        if(File.exist?(fasta_filename))
            @fasta_filename = fasta_filename
        else
            raise(ArgumentError,"ERROR: File '#{fasta_filename}' not found.")
        end
        if(min_len.class == Fixnum && min_len > 0)
            @min_len = min_len
        else
            raise(ArgumentError,"ERROR: '#{min_len}' must be a positive whole number.")
        end
    end
    def fasta_filename
        return @fasta_filename
    end
    def fasta_filehandl
        return @fasta_filehandl
    end
    def min_len
        return @min_len
    end
    def open()
        @fasta_filehandl = File.open(@fasta_filename,"r")
    end
    def nextSeq()
        seq_id = nil
        seq_bp_list = ""
        seq = nil
        if(!@fasta_filehandl.nil? && !@fasta_filehandl.closed? && !@fasta_filehandl.eof)
            begin
                line1 = @fasta_filehandl.gets
                if(line1.match(/^>(\w+)$/))
                    seq_id = $1
                    line2 = @fasta_filehandl.gets
                    if(line2.match(/^([atcgATCG]+)$/))
                        seq_bp_list = $1
                    end
                end
                if(!seq_id.nil? && !seq_bp_list.nil?)
                    if(seq_bp_list.length >= @min_len)
                        seq = Sequence.new(seq_id,seq_bp_list)
                    end
                end
            end while(seq.nil? && !@fasta_filehandl.eof)
        end
        return seq
    end
    def close()
        @fasta_filehandl.close
    end
end
