class NCBIBlaster
    CONNECTION_EXCEPTIONS = [Errno::ECONNRESET,Timeout::Error]
    NCBI_URI = URI('http://www.ncbi.nlm.nih.gov/blast/Blast.cgi')
    PutResponse = Struct.new(:rid,:seq)
    TEXT = Format.new("Text","txt")
    Format = Struct.new(:web_req_format,:file_suffix)
    def submitTblastxQuery(seq)
        return webCall(self.method(:put),seq)
    end
    def fetchTblastxResult(put_response)
        text_result = webCall(self.method(:get),TEXT,put_response)
        return buildNCBIResult(text_result,put_response.seq)
    end
    def buildNCBIResult(text_result,seq)
        ncbi_result = NCBIBlastResult.new(seq)
        if(!text_result.match(/No significant similarity found/))
            text_result.match(/\|([^\|]+)\|/)
            accession_num1 = $1
            accession_num2 = $2
            text_result.match(/\|[^\|]+\|.+[0-9][0-9]\.[0-9]\s+([0-9.-e]+)\s+/)
            e_value1 = $1
            e_value2 = $2
            align1 = Alignment.new(seq,accession_num1,e_value1)
            align2 = Alignment.new(seq,accession_num2,e_value2)
            ncbi_result.addAlignment(align1)
            ncbi_result.addAlignment(align2)
        end
        return ncbi_result
    end
    #pass method like self.method(:get)
    #pass params like TEXT,res
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
            put_result.body().match(/RID = ([0-9A-Z-]+)/)
            rid = $1
            if(rid)
                puts "RID: '#{rid}'"
                put_result.body().match(/RTOE = ([0-9]+)/)
                rtoe = $1
                puts "Estimated Request Execution Time: '#{rtoe}' seconds"
                return Put_Res.new(rid,seq.id,seq.bp_list)
            else
                puts "RID not returned. Retrying..."
                sleep(3)
            end
        end while(!put_result.body().match(/RID = [0-9A-Z-]+/))
    end
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
        end while(get_res_body.match(/Status=WAITING/))
        end_t = Time.now
        puts
        if(get_result.body().match(/Status=READY/))
            puts "Results completed after #{end_t - start_t} seconds"
            fn = "#{res.seq_name}.tblastx.#{format.file_suffix}"
            puts "Writing results to #{fn}..."
            fh = File.open(fn,"w")
            get_result.body.match(/(<p>.*?Query=)/m)
            res_header = $1
            res_content = get_result.body().gsub(res_header,'')
            fh.puts "Query ID=#{res.seq_name}\n#{res_content}"
            fh.close
            puts "done."
        else
            puts "UNKNOWN ERROR OCCURRED FOLLOWING GET"
        end
    end
end
