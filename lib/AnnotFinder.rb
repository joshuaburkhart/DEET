require 'net/http'
require_relative 'RgxLib'

class AnnotFinder
    NCBI_GENE_URI = URI('http://www.ncbi.nlm.nih.gov/gene')
    def initialize(acc_num,loghandl)
        if(!loghandl.nil? && loghandl.class == File && !loghandl.closed?)
            @loghandl = loghandl
            @acc_num = acc_num
            term_param = {
                :TERM => @acc_num
            }
            NCBI_GENE_URI.query = URI.encode_www_form(term_param)
            @annot_page = Net::HTTP.get_response(NCBI_GENE_URI)
        else
            raise(ArgumentError,"ERROR: loghandl '#{loghandl}' not a valid file handle")
        end
    end
    def getLocusTag
        if(!@annot_page.nil?)
            ans = nil
            @annot_page.body().match(RgxLib::ANNOT_LT_GRAB)
            locusTag = $1
            if(!locusTag.nil? && locusTag != "")
                ans = locusTag
            else
                msg = "ERROR: Locus tag for #{@acc_num} not found in annotation page\n#{@annot_page.body()}"
                @loghandl.puts msg
                ans = "EMPTY! (see log file '#{@loghandl.path}' for details)"
            end
            return ans
        else
            msg = "ERROR: annotation page not set"
            @loghandl.puts msg
            raise(StandardError,msg)
        end
    end
    def getName
        if(!@annot_page.nil?)
            ans = nil
            @annot_page.body().match(RgxLib::ANNOT_NM_GRAB)
            name = $1
            if(!name.nil? && name != "")
                ans = name
            else
                msg = "ERROR: name for #{@acc_num} not found in annotation page\n#{@annot_page.body()}"
                @loghandl.puts msg
                ans = "EMPTY! (see log file '#{@loghandl.path}' for details)"
            end
            return ans
        else
            msg = "ERROR: annotation page not set"
            @loghandl.puts msg
            raise(StandardError,msg)
        end
    end
end
