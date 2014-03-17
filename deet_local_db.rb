#!/usr/bin/ruby

require 'optparse'
require 'set'
require 'yaml'

require_relative 'lib/Sequence'
require_relative 'lib/Alignment'
require_relative 'lib/FastaParser'
require_relative 'lib/MicroArrayHashBuilder'
require_relative 'lib/RNAseqHashBuilder'
require_relative 'lib/RgxLib'
require_relative 'lib/LocalDbBlaster'
require_relative 'lib/LocalDbAnnotFinder'

START_TIME = Time.now
EST_OFFSET = (-5 * 3600)
EX_ID      = START_TIME.to_i
EST_TIME   = Time.now.localtime(EST_OFFSET)

Q_LIM        = 0.05 #for microarray
P_ADJ_LIM    = 0.05 #for rna seq
E_LIM        = 1.0e-5
MIN_SEQ_LEN  = 100

OUT_DIR_NAME = "output/"
SEQ_HEADER = "ID~Status~Name~Locus Tag~NCBI Acc #~Paralog #~Expr Sig"


%x(mkdir -p #{OUT_DIR_NAME})
outfile_prefix  = "#{OUT_DIR_NAME}#{EX_ID}"
result_filename = "#{outfile_prefix}.result"
resulthandl     = File.open(result_filename,"w")
seq_filename    = "#{outfile_prefix}.seqs.csv"
seqhandl        = File.open(seq_filename,"w")
log_filename    = "#{outfile_prefix}.log"
loghandl        = File.open(log_filename,"w")

msg = "INFO: Execution #{EX_ID} started.\nQ_LIM = #{Q_LIM}\nE_LIM = #{E_LIM}\nMIN_SEQ_LEN = #{MIN_SEQ_LEN}"
loghandl.puts msg
puts msg

options = {}
optparse = OptionParser.new { |opts|
    opts.banner = <<-EOS
Usage: ruby deet.rb -f <fasta file 1> ... <fasta file n> -m <dir containing only ma files>
    EOS
    opts.on('-h','--help','Display this screen'){
        puts opts
        exit
    }
    options[:fasta_files] = nil
    opts.on('-f','--fastas FASTA_FILE1,FASTA_FILE2',Array,'FASTA files') { |fasta_files|
        options[:fasta_files] = fasta_files
    }
    options[:ma_files] = nil
    opts.on('-m','--mad MA_DIR','Directory containing only valid MA files') { |mad|
        options[:ma_files] = Array.new
        Dir.entries(mad).each {|file|
            if(file != "." && file != "..")
                options[:ma_files] << "#{mad}/#{file}"
            end
        }
    }
    options[:rna_seq_files] = nil
    opts.on('-r','--rsqd RNA_SEQ_DIR','Directory containing only valid RNAseq files') { |rsqd|
        options[:rna_seq_files] = Array.new
        Dir.entries(rsqd).each {|file|
            if(file != "." && file != "..")
                options[:rna_seq_files] << "#{rsqd}/#{file}"
            end
        }
    }
}
optparse.parse!

def checkSeqsInHash(seqs,hash,hash_len,loghandl)
    msg = "INFO: Checking fasta sequences are in hashed data."
    loghandl.puts msg
    puts msg
    msg = "=================================================="
    loghandl.puts msg
    puts msg
    invalid_seqs = Set.new
    seqs.each {|seq|
        expr_sig = hash[seq.id]
        if(!expr_sig.nil? && expr_sig.class == String && expr_sig.match(RgxLib::ACCG_EXPR_SIG))
            print "."
            $stdout.flush
        else
            invalid_seqs << seq
            msg = "\nWARNING: sequence '#{seq.id}' not a valid key in hash."
            loghandl.puts msg
            puts msg
            $stdout.flush
            print "x"
            $stdout.flush
        end
    }
    puts
    seqs = seqs - invalid_seqs
    msg = "INFO: #{invalid_seqs.length} sequences removed, #{seqs.length} remain..."
    loghandl.puts msg
    puts msg
    return seqs
end

def parseBlastResultFromOutput(text_result,seq)
    ncbi_result = nil 
    if(!seq.nil? && !text_result.nil?)
        ncbi_result = NCBIBlastResult.new(seq)
        if(!text_result.match(/^#{seq.id}/))
            msg = "INFO: No match reported for #{seq.id}"
            puts msg
        else
            text_result.match(/^#{seq.id},(.*?),(.*?),/)
            accession_num = $1
            e_value = $2
            if(!accession_num.nil? && !e_value.nil?)
                align = Alignment.new(seq,accession_num,e_value)
                ncbi_result.addAlignment(align)
            else
                msg = "WARNING: Cannot parse alignment for #{seq.id}. accession_num=#{accession_num}, e_value=?#{e_value}"
                puts msg
            end
        end 
    else
        msg = "ERROR: NCBIBlastResult cannot be build using nil objects"
        raise(ArgumentError,msg)
    end 
    return ncbi_result
end


def parseFasta(fasta_files,id_grab_expr,loghandl,seqhandl)
    msg = "INFO: Below FASTA files supplied\n#{fasta_files.join("\n")}"
    loghandl.puts msg
    puts msg
    seqhandl.puts SEQ_HEADER
    seqs = Set.new
    fasta_files.each {|fasta_file|
        msg = "Loading #{fasta_file}..."
        loghandl.puts msg
        puts msg
        parser = FastaParser.new(fasta_file,id_grab_expr,loghandl)
        parser.open
        while(next_seq = parser.nextSeq)
            if(next_seq.bp_list.length > MIN_SEQ_LEN)
                seqs << next_seq
                print "."
            else
                seqhandl.puts next_seq
                print "x"
            end
            $stdout.flush
        end
        parser.close
    }
    puts
    msg = "INFO: Fasta files loaded. (#{seqs.length} unique sequences > #{MIN_SEQ_LEN}bp)"
    loghandl.puts msg
    puts msg
    msg = "=====================#{'='*seqs.length.to_s.length}================#{'='*MIN_SEQ_LEN.to_s.length}=" 
    loghandl.puts msg
    puts msg
    return seqs
end

if(options[:fasta_files].nil?)
    msg = "ERROR: No FASTA files supplied"
    loghandl.puts msg
    puts msg
    raise(ArgumentError,msg)
end

seq_hash = nil
fasta_seqs = nil
expr_sig_len = nil
if(options[:ma_files].nil? && options[:rna_seq_files].nil?)
    msg = "ERROR: No expression files supplied"
    loghandl.puts msg
    puts msg
    raise(ArgumentError,msg)
elsif(options[:rna_seq_files].nil?)
    msg = "INFO: Below MA files supplied (order maintained in expression signatures)\n#{options[:ma_files].join("\n")}"
    loghandl.puts msg
    puts msg
    expr_sig_len = options[:ma_files].length
    fasta_seqs = parseFasta(options[:fasta_files],RgxLib::FASTP_MA_ID_GRAB,loghandl,seqhandl)
    seq_hash = MicroArrayHashBuilder.makeHash(*options[:ma_files],Q_LIM,loghandl)
elsif(options[:ma_files].nil?)
    msg = "INFO: Below RNAseq files supplied (order maintained in expression signatures)\n#{options[:rna_seq_files].join("\n")}"
    loghandl.puts msg
    puts msg
    expr_sig_len = options[:rna_seq_files].length
    fasta_seqs = parseFasta(options[:fasta_files],RgxLib::FASTP_RNA_SEQ_ID_GRAB,loghandl,seqhandl)
    seq_hash = RNAseqHashBuilder.makeHash(*options[:rna_seq_files],P_ADJ_LIM,loghandl)
else
    msg = "ERROR: Only one type of expression file may be specified."
    loghandl.puts msg
    raise(ArgumentError,msg)
end
seqs = checkSeqsInHash(fasta_seqs,seq_hash,expr_sig_len,loghandl)

msg = "Querying Local Blast Database..."
loghandl.puts msg
puts msg
blaster            = LocalDbBlaster.new(loghandl)
local_db_blast_results = Hash.new
seq_keys           = Hash.new

puts "local_db_blast_results is a #{local_db_blast_results.class} with #{local_db_blast_results.count} members"
puts "seq_keys is a #{seq_keys.class} with #{seq_keys.count} members"

seq_batch = Set.new

QUERY_FILENAME = "query.fasta"
File.write(QUERY_FILENAME,"") #create an empty file

OUT_FILENAME = "out.txt"
if(File.exist?(OUT_FILENAME))
    File.delete(OUT_FILENAME) #assure output file does not exist
end

seq_batch_count = 0
seqs.each_with_index {|seq,i|
    File.write(QUERY_FILENAME,">#{seq.id}\n#{seq.bp_list}\n",File.size(QUERY_FILENAME),mode: 'a')
    seq_batch << seq
    if(i % 5000 == 4999)
        dt = Time.now - START_TIME
        h = (dt / 3600).floor
        m = ((dt % 3600) / 60).floor
        s = ((dt % 3600) % 60).floor
        printf("Submitting sequence batch #{seq_batch_count} to local blast db at T+%02.0f:%02.0f:%02.0f\n",h,m,s)        
        blaster.submitTblastxQuery(QUERY_FILENAME,OUT_FILENAME)
        while(!File.exist?(OUT_FILENAME))
            puts "Waiting for Sequence batch #{seq_batch_count} to finish blasting..."
            sleep(60)            
        end
        sleep(5) #wait for tblastx to stop writing output
        fh = File.open(OUT_FILENAME)
        text_result = fh.read
        fh.close
        seq_batch_seq_count = 1
        seq_batch.each{|seq_batch_seq|
            local_db_blast_result = parseBlastResultFromOutput(text_result,seq_batch_seq)
            if(!local_db_blast_result.nil?)
                local_db_blast_results[local_db_blast_result.sequence.id] = local_db_blast_result
                seq_keys[local_db_blast_result.sequence.id] = local_db_blast_result.sequence
                printf("Retrieved sequence #{(5000 * seq_batch_count) + seq_batch_seq_count}, #{local_db_blast_result.sequence.id}, from local blast db at T+%02.0f:%02.0f:%02.0f\n",h,m,s)
            else
                puts "Blast result nil!"
            end
            seq_batch_seq_count = seq_batch_seq_count + 1
        }
        seq_batch = Set.new
        File.write(QUERY_FILENAME,"") #clear query
        File.delete(OUT_FILENAME) #clear output file
        seq_batch_count = seq_batch_count + 1
    end
}
msg = "Query results retreived"
loghandl.puts msg
puts msg
msg = "======================="
loghandl.puts msg
puts msg

msg = "Grouping sequences..."
loghandl.puts msg
puts msg

acc_num_groups = Hash.new

local_db_blast_results.values.each_with_index {|local_db_res,j|
    print "."
    $stdout.flush
    if(local_db_res.hasAlignments?)
        acc_num = local_db_res.bestAlignment.accession_num        
        local_annot_finder = LocalDbAnnotFinder.new(acc_num,loghandl) #this searches the fasta used to create the local blast database
        annot_name = local_annot_finder.getName
        annot_locus_tag = local_annot_finder.getLocusTag
        acc_num_groups[acc_num] = AccessionNumGroup.new(annot_name,annot_locus_tag,acc_num,expr_sig_len,loghandl)
        seq_id   = local_db_res.sequence.id
        expr_sig = seq_hash[seq_id]
        if(!expr_sig.nil?)
            acc_num_groups[acc_num].addRes(local_db_res,expr_sig)
        else
            msg = "WARNING: expr_sig for #{seq_id} is nil. This may indicate FASTA/MA mismatch."
            loghandl.puts msg
            puts msg
        end
    else
        cur_seq          = local_db_res.sequence
        cur_seq.expr_sig = seq_hash[cur_seq.id]
        cur_seq.ignored  = false
        cur_seq.orphan   = true
        seqhandl.puts cur_seq.to_s
    end
}
puts
msg = "Sequences grouped by accession number and expression signature"
loghandl.puts msg
puts msg
msg = "=============================================================="
loghandl.puts msg
puts msg

msg = "Writing results..."
loghandl.puts msg
puts msg
acc_num_groups.values.each {|acc_num_group|
    print "."
    $stdout.flush
    seqhandl.puts(acc_num_group.repSeqDat(E_LIM))
    resulthandl.puts(acc_num_group.to_s)
    resulthandl.puts
}
puts
resulthandl.close
seqhandl.close
loghandl.close
puts "Result file: #{result_filename}"
puts "Seq file:    #{seq_filename}"
puts "Log file:    #{log_filename}"
