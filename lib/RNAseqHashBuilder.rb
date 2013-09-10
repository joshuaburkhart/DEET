require_relative 'RgxLib'

class RNAseqHashBuilder
    def self.makeHash(*rna_seq_filenames,p_adj_lim,loghandl)
        if(!loghandl.nil? && loghandl.class == File && !loghandl.closed?)
            @loghandl = loghandl
            msg = "Hashing RNAseq data..."
            @loghandl.puts msg 
            puts msg 
            if(!p_adj_lim.nil? && p_adj_lim.class == Float && p_adj_lim > 0.0)
                if(!rna_seq_filenames.nil? && rna_seq_filenames.class == Array && rna_seq_filenames.length > 0)
                    rna_seq_filenames.each {|f|
                        if(!File.exists?(f))
                            msg = "ERROR: file '#{f}' does not exist"
                            @loghandl.puts msg
                            raise(ArgumentError,msg)
                        end
                    }
                    seq_sigs = Hash.new
                    rna_seq_filenames.each {|f|
                        fhandl = File.open(f,"r")
                        header = fhandl.gets
                        if(!header.nil? && header.class == String && header.match(RgxLib::RSHB_HEADER_CHECK))
                            while(line = fhandl.gets)
                                if(!line.include?("NA"))
                                   line.match(RgxLib::RSHB_SEQ_ID_GRAB)
                                   seq_id = $1
                                   line.match(RgxLib::RSHB_P_ADJ_GRAB)
                                   p_adj = Float($1)
                                   sig_sample = "X" #not significant
                                   if(p_adj < p_adj_lim)
                                       line.match(RgxLib::RSHB_L2FC_GRAB)
                                       l2fc = $1 == "Inf" ? Float(1.0/0.0) : Float($1)
                                       sig_sample = l2fc > 0 ? "1" : "0"
                                   end
                                   if(seq_sigs[seq_id].nil?)
                                       seq_sigs[seq_id] = ""
                                   end
                                   seq_sigs[seq_id] << sig_sample
                                else
                                    msg = "WARNING: RNAseq file '#{f}' contains nonsense data line '#{line}'"
                                    @loghandl.puts msg
                                    puts msg
                                end
                            end
                        else
                            msg = "ERROR: file '#{f}' does not have a valid RNAseq header.\n#{header}"
                            @loghandl.puts msg
                            raise(ArgumentError,msg)
                        end
                    }
                    msg = "RNAseq data hashed"
                    @loghandl.puts msg 
                    puts msg 
                    msg = "======================"
                    @loghandl.puts msg 
                    puts msg 
                    return seq_sigs
                else
                    msg = "ERROR: rna_seq_filenames '#{rna_seq_filenames}' not a valid array of filenames"
                    @loghandl.puts msg
                    raise(ArgumentError,msg)
                end
            else
                msg = "ERROR: p_adj_lim '#{p_adj_lim}' not a valid adjusted p limit"
                @loghandl.puts msg
                raise(ArgumentError,msg)
            end
        else
            raise(ArgumentError,"ERROR: loghandl '#{loghandl}' not a valid file handle")
        end
    end
end
