require_relative 'RgxLib'

class MicroArrayHashBuilder
    def self.makeHash(*ma_filenames,q_lim,loghandl)
        if(!loghandl.nil? && loghandl.class == File && !loghandl.closed?)
            @loghandl = loghandl
            msg = "Hashing microarray data..."
            @loghandl.puts msg
            puts msg
            if(!q_lim.nil? && q_lim.class == Float && q_lim > 0.0)
                if(!ma_filenames.nil? && ma_filenames.class == Array && ma_filenames.length > 0)
                    ma_filenames.each {|f|
                        if(!File.exists?(f))
                            msg = "ERROR: file '#{f}' does not exist"
                            @loghandl.puts msg
                            raise(ArgumentError,msg)
                        end
                    }
                    seq_sigs = Hash.new
                    ma_filenames.each {|f|
                        fhandl = File.open(f,"r")
                        header = fhandl.gets
                        if(!header.nil? && header.class == String && header.match(RgxLib::MAHB_HEADER_CHECK))
                            while(line = fhandl.gets)
                                line.match(RgxLib::MAHB_SEQ_ID_GRAB)
                                seq_id = $1
                                line.match(RgxLib::MAHB_Q_VAL_GRAB)
                                q_val = Float($9)
                                sig_sample = "X" #not significant
                                if(q_val < q_lim)
                                    line.match(RgxLib::MAHB_M_VAL_GRAB)
                                    m_val = Float($1)
                                    sig_sample = m_val > 0 ? "1" : "0"
                                end
                                if(seq_sigs[seq_id].nil?)
                                    seq_sigs[seq_id] = ""
                                end
                                seq_sigs[seq_id] << sig_sample
                            end
                        else
                            msg = "ERROR: file '#{f}' does not have a valid MA header.\n#{header}"
                            @loghandl.puts msg
                            raise(ArgumentError,msg)
                        end
                    }
                    msg = "Microarray data hashed"
                    @loghandl.puts msg 
                    puts msg 
                    msg = "======================"
                    @loghandl.puts msg 
                    puts msg
                    return seq_sigs
                else
                    msg = "ERROR: ma_filenames '#{ma_filenames}' not a valid array of filenames"
                    @loghandl.puts msg
                    raise(ArgumentError,msg)
                end
            else
                msg = "ERROR: q_lim '#{q_lim}' not a valid q value limit"
                @loghandl.puts msg
                raise(ArgumentError,msg)
            end
        else
            raise(ArgumentError,"ERROR: loghandl '#{loghandl}' not a valid file handle")
        end
    end
end
