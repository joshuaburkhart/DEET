require_relative 'RgxLib'

class FastaParser
    @fasta_filename
    @fasta_filehandl
    def initialize(fasta_filename,loghandl)
        if(!loghandl.nil? && loghandl.class == File && !loghandl.closed?)
            @loghandl = loghandl
            if(!fasta_filename.nil? && fasta_filename.class == String && File.exist?(fasta_filename))
                @fasta_filename = fasta_filename
            else
                msg = "ERROR: File '#{fasta_filename}' not a valid file"
                @loghandl.puts msg
                raise(ArgumentError,msg)
            end
        else
            raise(ArgumentError,"ERROR: loghandl '#{loghandl}' not a valid file handle")
        end
    end
    def fasta_filename
        return @fasta_filename
    end
    def fasta_filehandl
        return @fasta_filehandl
    end
    def open()
        if(!@fasta_filename.nil?)
            @fasta_filehandl = File.open(@fasta_filename,"r")
        else
            msg = "File not specified"
            @loghandl.puts msg
            raise(IOError,msg)
        end
    end
    def nextSeq()
        seq_id = nil
        seq_bp_list = ""
        seq = nil
        id_line = nil
        while(seq.nil? && !fileHandlIssue)
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
                    seq = Sequence.new(seq_id,seq_bp_list)
                end
            end
        end
        return seq
    end
    def readFastaLine
        line = nil
        if(!fileHandlIssue)
            line = @fasta_filehandl.gets
            if(!line.match(RgxLib::FASTP_ID_GRAB) && !line.match(RgxLib::SEQ_BP_LIST))
                msg = "ERROR: '#{line}' not valid fasta format"
                @loghandl.puts msg
                raise(ArgumentError,msg)
            end
        end
        return line
    end
    def fileHandlIssue
        return @fasta_filehandl.nil? || @fasta_filehandl.closed? || @fasta_filehandl.eof
    end
    def close()
        @fasta_filehandl.close
    end
end
