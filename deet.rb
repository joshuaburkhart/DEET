#!/usr/bin/ruby

require 'optparse'
require 'set'
require 'yaml'

require_relative 'lib/Sequence'
require_relative 'lib/Alignment'
require_relative 'lib/NCBIBlastResult'
require_relative 'lib/NCBIBlaster'
require_relative 'lib/AccessionNumGroup'
require_relative 'lib/FastaParser'
require_relative 'lib/MicroArrayHashBuilder'
require_relative 'lib/RNAseqHashBuilder'
require_relative 'lib/RgxLib'

START_TIME = Time.now
EST_OFFSET = (-5 * 3600)
EX_ID      = START_TIME.to_i
EST_TIME   = Time.now.localtime(EST_OFFSET)

Q_LIM        = 0.05 #for microarray
P_ADJ_LIM    = 0.05 #for rna seq
E_LIM        = 1.0e-5
MIN_SEQ_LEN  = 100

OUT_DIR_NAME = "output/"
NCBI_RES_SER_FN = "ncbi_blast_results.serialized"
SEQ_KEYS_SER_FN = "seq_keys.serialized"
ACC_NUMG_SER_FN = "acc_num_groups.serialized"
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

if(!EST_TIME.saturday? && !EST_TIME.sunday? && !(EST_TIME.hour > 21))
    puts "Due to the intensive load this program may put on NCBI servers, this program should only be run on weekends or between the hours of 9pm and 5am."
    puts "Choose One:"
    puts "w - wait for next eligible execution time, then run automatically"
    puts "a - abort execution and manually run this program later, or not"
    puts "r - run immediately, ignoring NCBI guidelines"
    print ">"
    answer = $stdin.gets.strip
    if(!answer.nil? && answer.class == String && answer.length == 1 && answer.match(/[war]/))
        if(answer == "w")
            num_minutes = 0
            t = Time.now.localtime(EST_OFFSET)
            puts "Waiting..."
            while(!t.saturday? && !t.sunday? && !(t.hour > 21))
                print "."
                $stdout.flush
                sleep(600) #10 minutes
                t = Time.now.localtime(EST_OFFSET)
            end         
            msg = "Beginning execution after #{num_minutes} minute wait..."
            loghandl.puts msg
            puts msg
            run = true
        elsif(answer == "a")
            puts "Aborted."
            exit(0)
        elsif(answer == "r")
            msg = "Running at discouraged time. Good luck..."
            loghandl.puts msg
            puts msg
        end
    end
end

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

def storePutResult(ret_seq_count,p_res,blaster,ncbi_blast_results,seq_keys)
    dt = Time.now - START_TIME
    h  = (dt / 3600).floor
    m  = ((dt % 3600) / 60).floor
    s  = ((dt % 3600) % 60).floor
    printf("Retrieving sequence #{ret_seq_count}, #{p_res.seq.id}, from NCBI at T+%02.0f:%02.0f:%02.0f\n",h,m,s)
    ncbi_blast_result = blaster.fetchTblastxResult(p_res)
    if(!ncbi_blast_result.nil?)
        ncbi_blast_results[ncbi_blast_result.sequence.id] = ncbi_blast_result
        seq_keys[ncbi_blast_result.sequence.id] = ncbi_blast_result.sequence
        puts "Serializing #{ncbi_blast_results.size} blast results..."
        serialized_results    = YAML::dump(ncbi_blast_results)
        serialized_keys       = YAML::dump(seq_keys)
        serialized_results_fh = File.open(NCBI_RES_SER_FN,"w")
        serialized_keys_fh    = File.open(SEQ_KEYS_SER_FN,"w")
        serialized_results_fh.write(serialized_results)
        serialized_keys_fh.write(serialized_keys)
        serialized_results_fh.flush
        serialized_keys_fh.flush
        serialized_results_fh.close
        serialized_keys_fh.close
        puts "Blast results serialized..."
    else
        puts "Blast result nil!"
    end
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

msg = "Querying NCBI..."
loghandl.puts msg
puts msg
blaster            = NCBIBlaster.new(loghandl)
put_results        = Set.new
ncbi_blast_results = nil
seq_keys = nil
if(File.exists?(NCBI_RES_SER_FN) && File.exists?(SEQ_KEYS_SER_FN))
    puts "Detected serialized blast results..."
    serialized_results = File.read(NCBI_RES_SER_FN)
    serialized_keys    = File.read(SEQ_KEYS_SER_FN)
    ncbi_blast_results = YAML::load(serialized_results)
    seq_keys           = YAML::load(serialized_keys)
    puts "Loaded #{ncbi_blast_results.size} blast results..."
end
if(ncbi_blast_results.nil? || seq_keys.nil?)
    ncbi_blast_results = Hash.new
    seq_keys           = Hash.new
end
puts "ncbi_blast_results is a #{ncbi_blast_results.class} with #{ncbi_blast_results.count} members"
puts "seq_keys is a #{seq_keys.class} with #{seq_keys.count} members"

ret_seq_count      = 0
seqs.each_with_index {|seq,i|
    if(!seq_keys.key?(seq.id))
        dt = Time.now - START_TIME
        h = (dt / 3600).floor
        m = ((dt % 3600) / 60).floor
        s = ((dt % 3600) % 60).floor
        printf("Submitting sequence #{i}, id=#{seq.id}, bp_list=#{seq.bp_list}, to NCBI at T+%02.0f:%02.0f:%02.0f\n",h,m,s)
        result = blaster.submitTblastxQuery(seq)
        if(!result.nil?)
            put_results << result
        end
        if(i % 100 == 99)
            put_results.each {|p_res|
                storePutResult(ret_seq_count,p_res,blaster,ncbi_blast_results,seq_keys)
                put_results.delete(p_res)
                ret_seq_count += 1
            }
        end
    else
        puts "seq #{seq} already in blast results!!!"
    end
}
if(put_results.length > 0)
    put_results.each {|rem_p_res|
        storePutResult(ret_seq_count,rem_p_res,blaster,ncbi_blast_results,seq_keys)
        put_results.delete(rem_p_res)
        ret_seq_count += 1
    }
end
msg = "Query results retreived"
loghandl.puts msg
puts msg
msg = "======================="
loghandl.puts msg
puts msg

msg = "Grouping sequences..."
loghandl.puts msg
puts msg
acc_num_groups = nil
if(File.exists?(ACC_NUMG_SER_FN))
    puts "Detected serialized accession number groups..."
    serialized_acc_num_groups = File.read(ACC_NUMG_SER_FN)
    acc_num_groups = YAML::load(serialized_acc_num_groups)
    puts "Loaded #{acc_num_groups.size} groups..."
end
if(acc_num_groups.nil?)
    acc_num_groups = Hash.new
end
ncbi_blast_results.values.each_with_index {|ncbi_res,j|
    print "."
    $stdout.flush
    if(ncbi_res.hasAlignments?)
        acc_num = ncbi_res.bestAlignment.accession_num        
        if(acc_num_groups[acc_num].nil?)
            acc_num_groups[acc_num] = AccessionNumGroup.new(acc_num,expr_sig_len,loghandl)
        else
            puts "acc_num #{acc_num} already in acc_num_groups!!!"
        end
        seq_id   = ncbi_res.sequence.id
        expr_sig = seq_hash[seq_id]
        if(!expr_sig.nil?)
            acc_num_groups[acc_num].addRes(ncbi_res,expr_sig)
            puts "Serializing #{acc_num_groups.size} accession number groups..."
            serialized_acc_num_groups = YAML::dump(acc_num_groups)
            serialized_acc_num_groups_fh = File.open(ACC_NUMG_SER_FN,"w")
            serialized_acc_num_groups_fh.write(serialized_acc_num_groups)
            serialized_acc_num_groups_fh.flush
            serialized_acc_num_groups_fh.close
            puts "Accession number groups serialized..."
        else
            msg = "WARNING: expr_sig for #{seq_id} is nil. This may indicate FASTA/MA mismatch."
            @loghandl.puts msg
            puts msg
        end
    else
        cur_seq          = ncbi_res.sequence
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
