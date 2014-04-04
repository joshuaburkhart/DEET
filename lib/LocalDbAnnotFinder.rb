require 'net/http'
require_relative 'RgxLib'
require_relative 'NCBIBlaster'

class LocalDbAnnotFinder < NCBIBlaster
    LOCAL_DB_FASTA_FILENAME = "/N/u/joshburk/Mason/refseq_complete/invertebrate_rna.fna.headers_only"
    @name
    @locus_tag
    def initialize(acc_num,loghandl)
        if(!loghandl.nil? && loghandl.class == File && !loghandl.closed?)
            @loghandl = loghandl
            @acc_num = acc_num
            fetchAnnotation
        else
            raise(ArgumentError,"ERROR: loghandl '#{loghandl}' not a valid file handle")
        end
        super(loghandl)
    end
    def fetchAnnotation
        @annot_line = %x(cat #{LOCAL_DB_FASTA_FILENAME} | grep #{@acc_num})
    end
    def getLocusTag
        if(@locus_tag.nil?)
            if(!@annot_line.nil?)
                ans = nil
                @annot_line.match(/\|.*\|.*\|\s*.*\s*\((.*)\).*/)
                locusTag = $1
                if(!locusTag.nil? && locusTag != "")
                    ans = locusTag
                else
                    msg = "ERROR: Locus tag for #{@acc_num} not found in annotation line\n#{@annot_line}"
                    @loghandl.puts msg
                    ans = "EMPTY! (see log file '#{@loghandl.path}' for details)"
                end
                @locus_tag = ans
            else
                msg = "ERROR: annotation not set"
                @loghandl.puts msg
                raise(StandardError,msg)
            end
        end
        return @locus_tag
    end
    def getName
        if(@name.nil?)
            if(!@annot_line.nil?)
                ans = nil
                @annot_line.match(/\|.*\|.*\|\s*(.*)\s*\.*\.*/)
                name = $1
                if(!name.nil? && name != "")
                    ans = name
                else
                    msg = "ERROR: name for #{@acc_num} not found in annotation line\n#{@annot_line}"
                    @loghandl.puts msg
                    ans = "EMPTY! (see log file '#{@loghandl.path}' for details)"
                end
                @name = ans
            else
                msg = "ERROR: annotation line nil for acc num #{@acc_num}"
                @loghandl.puts msg
                raise(StandardError,msg)
            end
        end
        return @name
    end
end
