require 'net/http'
require_relative 'RgxLib'
require_relative 'NCBIBlaster'

class AnnotFinder < NCBIBlaster
    NCBI_GENE_URI = URI('http://www.ncbi.nlm.nih.gov/gene')
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
        @annot_page = webCall(self.method(:annotate))
    end
    def annotate
        term_param = {
            :TERM => @acc_num
        }
        NCBI_GENE_URI.query = URI.encode_www_form(term_param)
        ncbi_response = Net::HTTP.get_response(NCBI_GENE_URI)
        return ncbi_response.nil? ? "<ERROR QUERYING NCBI>" : ncbi_response
    end
    def getLocusTag
        if(@locus_tag.nil?)
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
                @locus_tag = ans
            else
                msg = "ERROR: annotation page not set"
                @loghandl.puts msg
                raise(StandardError,msg)
            end
        end
        return @locus_tag
    end
    def getName
        if(@name.nil?)
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
                @name = ans
            else
                msg = "ERROR: annotation page not set"
                @loghandl.puts msg
                raise(StandardError,msg)
            end
        end
        return @name
    end
end
