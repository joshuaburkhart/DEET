require 'net/http'
require_relative 'Sequence'
require_relative 'Alignment'
require_relative 'NCBIBlastResult'

class NCBIBlaster
    CONNECTION_EXCEPTIONS = [Errno::ECONNRESET,Timeout::Error,EOFError]
    NCBI_URI = URI('http://www.ncbi.nlm.nih.gov/blast/Blast.cgi')
    PutResponse = Struct.new(:rid,:seq)
    Format = Struct.new(:web_req_format,:file_suffix)
    TEXT = Format.new("Text","txt")
    T_LIM = 10 * 60
    def initialize(loghandl)
        if(!loghandl.nil? && loghandl.class == File && !loghandl.closed?)
            @loghandl = loghandl
            msg = "INFO: NCBIBlaster initialized as '#{self.to_s}'"
            @loghandl.puts msg
        else
            raise(ArgumentError,"ERROR: loghandl '#{loghandl}' not a valid file handle")
        end
    end
    def submitTblastxQuery(seq)
        return webCall(self.method(:put),seq)
    end
    def fetchTblastxResult(put_response)
        blast_result = nil
        if(!put_response.nil?)
            if(!put_response.seq.nil?)
                text_result = webCall(self.method(:get),TEXT,put_response)
                if(!text_result.nil?)
                    blast_result = buildNCBIResult(text_result,put_response.seq)
                end
            end
        else
            msg = "ERROR: put_response nil"
            @loghandl.puts msg
            raise(ArgumentError,msg)
        end
        return blast_result
    end
    def buildNCBIResult(text_result,seq)
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
                    text_result.match(RgxLib::BLST_E_VAL_GRAB)
                    e_value = $1
                    if(!accession_num.nil? && !e_value.nil?)
                        align = Alignment.new(seq,accession_num,e_value)
                        ncbi_result.addAlignment(align)
                    else
                        msg = "WARNING: Cannot parse alignment for #{seq.id}\n#{text_result}" 
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
        cur_t = Time.now
        begin
            get_result = Net::HTTP.get_response(NCBI_URI)
            print "."
            STDOUT.flush
            get_res_body = get_result.body()
            sleep(3)
            cur_t = Time.now
        end while(get_res_body.match(RgxLib::BLST_WAIT) && (cur_t - start_t < T_LIM))
        end_t = Time.now
        puts
        if(get_result.body().match(RgxLib::BLST_READY)
           puts "Results completed after #{end_t - start_t} seconds"
           get_result.body.match(RgxLib::BLST_HEADER_GRAB))
           res_header = $1
           res_content = get_result.body().gsub(res_header,'')
           res_text = "Query ID=#{res.seq.id}\n#{res_content}"
           return res_text
        elsif(cur_t - start_t >= T_LIM)
            msg = "INFO: get request for #{res} exceeded time limit (#{T_LIM} seconds)"
            @loghandl.puts msg
        else
            msg = "WARNING: unrecognized response for #{res} following get request\n#{get_result.body()}"
            @loghandl.puts msg
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
                    sleep(60)
                    retry
                else
                    msg = "WARNING: Maximum retry attempts reached, unable to get results for #{params.join(', ')}"
                    @loghandl.puts msg
                end 
            end 
        else
            msg = "wrong number of arguments (#{args_supplied} for #{args_required})"
            @loghandl.puts msg
            raise(ArgumentError,msg)
        end 
    end
end
