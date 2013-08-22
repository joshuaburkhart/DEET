require_relative 'RgxLib'

class MicroArrayHashBuilder
    QLIM = 0.05
    def self.makeHash(*ma_filenames)
        if(!ma_filenames.nil? && ma_filenames.class == Array && ma_filenames.length > 0)
            ma_filenames.each {|f|
                if(!File.exists?(f))
                    raise(ArgumentError,"ERROR: file '#{f}' does not exist")
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
                    %x(cat #{f})
                    raise(ArgumentError,"ERROR: file '#{f}' not valid ma format")
                end
            }
            return seq_sigs
        else
            raise(ArgumentError,"ERROR: ma_filenames '#{ma_filenames}' not a valid array of filenames")
        end
    end
end
