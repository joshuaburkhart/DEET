require_relative 'RgxLib'

class SimpleFastaParser
    @fasta_filename
    @fasta_filehandl
    @min_len
    def initialize(fasta_filename,min_len)
        if(!fasta_filename.nil? && fasta_filename.class == String && File.exist?(fasta_filename))
            if(!min_len.nil? && min_len.class == Fixnum && min_len > 0)
                @fasta_filename = fasta_filename
                @min_len = min_len
            else
                raise(ArgumentError,"ERROR: min_len '#{min_len}' must be a positive whole number")
            end
        else
            raise(ArgumentError,"ERROR: File '#{fasta_filename}' not a valid file")
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
        if(!@fasta_filename.nil?)
            @fasta_filehandl = File.open(@fasta_filename,"r")
        else
            raise(IOError,"File not specified")
        end
    end
    def nextSeq()
        seq_id = nil
        seq_bp_list = ""
        seq = nil
        id_line = nil
        if(@buffer.nil?)
            id_line = readFastaLine
        else
            id_line = @buffer
        end
        if(!id_line.nil?)
            if(id_line.match(RgxLib::FASTP_ID_GRAB))
                seq_id = $1
                bp_line = readFastaLine
                while(!bp_line.nil? && bp_line.match(RgxLib::SEQ_BP_LIST))
                    seq_bp_list += $1
                    bp_line = readFastaLine
                end
                @buffer = bp_line
            end
            if(!seq_id.nil? && !seq_bp_list.nil?)
                if(seq_bp_list.length >= @min_len)
                    seq = Sequence.new(seq_id,seq_bp_list)
                end
            end
        end
        return seq
    end
    def readFastaLine
        line = nil
        if(!@fasta_filehandl.nil? && !@fasta_filehandl.closed? && !@fasta_filehandl.eof)
            line = @fasta_filehandl.gets
            if(!line.match(RgxLib::FASTP_ID_GRAB) && !line.match(RgxLib::SEQ_BP_LIST))
                raise(ArgumentError,"ERROR: '#{line}' not valid fasta format")
            end
        end
        return line
    end
    def close()
        @fasta_filehandl.close
    end
end
