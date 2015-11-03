
require_relative 'Sequence'
require_relative 'Alignment'
require_relative 'NCBIBlastResult'

class LocalDbBlaster

    def initialize(loghandl)
        if(!loghandl.nil? && loghandl.class == File && !loghandl.closed?)
            @loghandl = loghandl
            msg = "INFO: LocalDbBlaster initialized as '#{self.to_s}'"
            @loghandl.puts msg
        else
            raise(ArgumentError,"ERROR: loghandl '#{loghandl}' not a valid file handle")
        end
    end
    def submitTblastxQuery(query_filename,out_filename)
        %x(/N/u/joshburk/Mason/blast/ncbi-blast-2.2.29+-src/c++/ReleaseMT/bin/tblastx -db /N/u/joshburk/Mason/refseq_complete/invertebrate_rna.fna -evalue 1.0e-5 -query #{query_filename} -num_threads 32 -out #{out_filename} -outfmt "10 qseqid sacc evalue" -max_target_seqs 2 -culling_limit 1)
    end
    def buildLocalDbBlastResult(text_result,seq)
        ncbi_result = nil
        if(!seq.nil? && !text_result.nil?)
            ncbi_result = NCBIBlastResult.new(seq)
            if(text_result.match(RgxLib::BLST_NO_MATCH))
                ncbi_result.valid = false
                msg = "INFO: No match reported for #{seq.id}\n#{text_result}"
                @loghandl.puts msg
            else
                2.times {
                    text_result.match(RgxLib::BLST_ACCN_GRAB)
                    accession_num = $1
                    text_result.match(/\|[^\|]+\|.+\<\/a\>\s+(\S+)/)
                    e_value = $1
                    if(!accession_num.nil? && !e_value.nil?)
                        align = Alignment.new(seq,accession_num,e_value)
                        ncbi_result.addAlignment(align)
                    else
                        msg = "WARNING: Cannot parse alignment for #{seq.id}. accession_num=#{accession_num}, e_value=?#{e_value}\n#{text_result}"
                        @loghandl.puts msg
                    end
                }
            end
        else
            msg = "ERROR: NCBIBlastResult cannot be build using nil objects"
            @loghandl.puts msg
            raise(ArgumentError,msg)
        end
        return ncbi_result
    end
end
