require_relative 'RgxLib'

class MicroArrayHashBuilder
    QLIM = 0.05
    def self.makeHash(*ma_filenames,loghandl)
        if(!loghandl.nil? && loghandl.class == File && !loghandl.closed?)
            if(!ma_filenames.nil? && ma_filenames.class == Array && ma_filenames.length > 0)
                ma_filenames.each {|f|
                    if(!File.exists?(f))
                        msg = "ERROR: file '#{f}' does not exist"
                        loghandl.puts msg
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
                            if(q_val < QLIM)
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
                        msg = "ERROR: file '#{f}' do not have a valid MA header"
                        loghandl.puts msg
                        raise(ArgumentError,msg)
                    end
                }
                return seq_sigs
            else
                msg = "ERROR: ma_filenames '#{ma_filenames}' not a valid array of filenames"
                loghandl.puts msg
                raise(ArgumentError,msg)
            end
        else
            raise(ArgumentError,"ERROR: loghandl '#{loghandl}' not a valid file handle")
        end
    end
end
