require 'net/http'
require_relative 'Sequence'
require_relative 'Alignment'
require_relative 'NCBIBlastResult'

class NCBIBlaster
    CONNECTION_EXCEPTIONS = [Errno::ECONNRESET,Timeout::Error]
    NCBI_URI = URI('http://www.ncbi.nlm.nih.gov/blast/Blast.cgi')
    PutResponse = Struct.new(:rid,:seq)
    Format = Struct.new(:web_req_format,:file_suffix)
    TEXT = Format.new("Text","txt")
    def submitTblastxQuery(seq)
        return webCall(self.method(:put),seq)
    end
    def fetchTblastxResult(put_response)
        if(!put_response.nil?)
            text_result = webCall(self.method(:get),TEXT,put_response)
            blast_result = buildNCBIResult(text_result,put_response.seq)
        else
            raise(ArgumentError,"ERROR: put_response nil")
        end
        return blast_result
    end
    def buildNCBIResult(text_result,seq)
        if(seq && text_result)
            ncbi_result = NCBIBlastResult.new(seq)
            if(text_result.match(RgxLib::BLST_NO_MATCH))
                ncbi_result.valid = false
            else
                2.times {
                    text_result.match(RgxLib::BLST_ACCN_GRAB)
                    accession_num = $1
                    text_result.match(RgxLib::BLST_E_VAL_GRAB)
                    e_value = $1
                    align = Alignment.new(seq,accession_num,e_value)
                    ncbi_result.addAlignment(align)
                }
            end
        else
            raise(ArgumentError,"ERROR: NCBIBlastResult cannot be build using nil objects")
        end
        return ncbi_result
    end
    private
    def put(seq)
        put_params = {
            :QUERY => seq.bp_list,
            :DATABASE => "nr",
            :HITLIST_SIZE => 10,
            :FILTER => 'L',
            :EXPECT => 10,
            :FORMAT_TYPE => "HTML",
            :PROGRAM => "tblastx",
            :CLIENT => "web",
            :SERVICE => "plain",
            :NCBI_GI => "on",
            :PAGE => "Nucleotides",
            :CMD => "Put",
        }
        put_result = nil
        begin
            NCBI_URI.query = URI.encode_www_form(put_params)
            put_result = Net::HTTP.get_response(NCBI_URI)
            put_result.body().match(RgxLib::BLST_RID_GRAB)
            rid = $1
            if(rid)
                put_result.body().match(RgxLib::BLST_RTOE_GRAB)
                rtoe = $1
                puts "Estimated Request Execution Time: '#{rtoe}' seconds"
                return PutResponse.new(rid,seq)
            else
                puts "RID not returned. Retrying..."
                sleep(60)
            end
        end while(!put_result.body().match(RgxLib::BLST_RID_GRAB))
    end
    private
    def get(format,res)
        get_params = {
            :RID => res.rid,
            :FORMAT_OBJECT => "Alignment",
            :FORMAT_TYPE => format.web_req_format,
            :DESCRIPTIONS => 10,
            :ALIGNMENTS => 10,
            :ALIGNMENT_TYPE => "Pairwise",
            :OVERVIEW => "yes",
            :CMD => "Get",
        }
        NCBI_URI.query = URI.encode_www_form(get_params)
        get_res_body = nil
        start_t = Time.now
        begin
            get_result = Net::HTTP.get_response(NCBI_URI)
            print "."
            STDOUT.flush
            get_res_body = get_result.body()
            sleep(1)
        end while(get_res_body.match(RgxLib::BLST_WAIT))
        end_t = Time.now
        puts
        if(get_result.body().match(RgxLib::BLST_READY)
           puts "Results completed after #{end_t - start_t} seconds"
           get_result.body.match(RgxLib::BLST_HEADER_GRAB))
           res_header = $1
           res_content = get_result.body().gsub(res_header,'')
           res_text = "Query ID=#{res.seq.id}\n#{res_content}"
           return res_text
        else
            puts "UNKNOWN ERROR OCCURRED FOLLOWING GET"
        end
    end
    #pass method like self.method(:get)
    #pass params like TEXT,res
    private
    def webCall(method,*params)
        args_supplied = params.length
        args_required = method.arity
        if(args_required == args_supplied)
            attempts = 0 
            begin
                method[*params]
            rescue *CONNECTION_EXCEPTIONS
                attempts += 1
                puts "Recovered from connection error..."
                if(attempts < 10) 
                    puts "Waiting to retry..."
                    sleep(30)
                    retry
                else
                    puts "Unable to get results for #{params.join(', ')}."
                end 
            end 
        else
            raise ArgumentError, "wrong number of arguments (#{args_supplied} for #{args_required})"
        end 
    end
end
