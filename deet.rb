#!/usr/bin/ruby

require 'optparse'
require 'set'

require_relative 'lib/Sequence'
require_relative 'lib/Alignment'
require_relative 'lib/NCBIBlastResult'
require_relative 'lib/NCBIBlaster'
require_relative 'lib/AccessionNumGroup'
require_relative 'lib/FastaParser'
require_relative 'lib/MicroArrayHashBuilder'
require_relative 'lib/RgxLib'

START_TIME = Time.now
EST_OFFSET = (-5 * 3600)
EX_ID      = START_TIME.to_i
EST_TIME   = Time.now.localtime(EST_OFFSET)

Q_LIM        = 0.05
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
Usage: ruby deet.rb -f <fasta file 1> ... <fasta file n> -m <ma file 1> ... <ma file n>

Example: ruby deet.rb -f input_files/fastas/singletons.self.remained.fna -m input_files/mas/photoperiod/limma.KCLD10-KCLD22.gene.de.txt,input_files/mas/photoperiod/limma.KCLD22-KCSD22.gene.de.txt,input_files/mas/photoperiod/limma.KCSD10-KCLD10.gene.de.txt,input_files/mas/photoperiod/limma.KCSD22-KCSD10.gene.de.txt,input_files/mas/photoperiod/limma.PBLD10-KCLD10.gene.de.txt,input_files/mas/photoperiod/limma.PBLD10-PBLD22.gene.de.txt,input_files/mas/photoperiod/limma.PBLD22-KCLD22.gene.de.txt,input_files/mas/photoperiod/limma.PBSD10-KCSD10.gene.de.txt,input_files/mas/photoperiod/limma.PBSD10-PBLD10.gene.de.txt,input_files/mas/photoperiod/limma.PBSD22-KCSD22.gene.de.txt,input_files/mas/photoperiod/limma.PBSD22-PBLD22.gene.de.txt,input_files/mas/photoperiod/limma.PBSD22-PBSD10.gene.de.txt
    EOS
    opts.on('-h','--help','Display this screen'){
        puts opts
        exit
    }
    options[:fasta_files] = nil
    opts.on('-f','--fastas FASTA_FILE1,FASTA_FILE2',Array,'FASTA files') { |fasta_files|
        options[:fasta_files] = fasta_files
        msg = "Accepted -f '#{options[:fasta_files]}' as argument."
        loghandl.puts msg
        puts msg
    }
    options[:ma_files] = nil
    opts.on('-m','--mas MA_FILE1,MA_FILE2',Array,'MA files') { |ma_files|
        options[:ma_files] = ma_files
        msg = "Accepted -m '#{options[:ma_files]}' as argument."
        loghandl.puts msg
        puts msg
    }
}
optparse.parse!

if(options[:fasta_files].nil?)
    msg = "ERROR: No FASTA files supplied"
    loghandl.puts msg
    puts msg
    raise(ArgumentError,msg)
else
    msg = "INFO: Below FASTA files supplied\n#{options[:fasta_files].join("\n")}"
    loghandl.puts msg
    puts msg
end
if(options[:ma_files].nil?)
    msg = "ERROR: No MA files supplied"
    loghandl.puts msg
    puts msg
    raise(ArgumentError,msg)
else
    msg = "INFO: Below MA files supplied (order maintained in expression signatures)\n#{options[:ma_files].join("\n")}"
    loghandl.puts msg
    puts msg
end

msg = "Hashing microarray data..."
loghandl.puts msg
puts msg
seq_hash = MicroArrayHashBuilder.makeHash(*options[:ma_files],Q_LIM,loghandl)
msg = "Microarray data hashed"
loghandl.puts msg
puts msg
msg = "======================"
loghandl.puts msg
puts msg

seqhandl.puts SEQ_HEADER
seqs = Set.new
options[:fasta_files].each {|fasta_file|
    msg = "Loading #{fasta_file}..."
    loghandl.puts msg
    puts msg
    parser = FastaParser.new(fasta_file,loghandl)
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
msg = "Fasta files loaded. (#{seqs.length} unique sequences > #{MIN_SEQ_LEN}bp)"
loghandl.puts msg
puts msg
msg = "=====================#{'='*seqs.length.to_s.length}================#{'='*MIN_SEQ_LEN.to_s.length}=" 
loghandl.puts msg
puts msg

msg = "Querying NCBI..."
loghandl.puts msg
puts msg
blaster            = NCBIBlaster.new(loghandl)
put_results        = Set.new
ncbi_blast_results = Set.new
ret_seq_count      = 0
seqs.each_with_index {|seq,i|
    dt = Time.now - START_TIME
    h = (dt / 3600).floor
    m = ((dt % 3600) / 60).floor
    s = ((dt % 3600) % 60).floor
    printf("Submitting sequence #{i}, #{seq.id}, to NCBI at T+%02.0f:%02.0f:%02.0f\n",h,m,s)
    put_results << blaster.submitTblastxQuery(seq)
    if(i % 100 == 99)
        put_results.each {|p_res|
            dt = Time.now - START_TIME
            h  = (dt / 3600).floor
            m  = ((dt % 3600) / 60).floor
            s  = ((dt % 3600) % 60).floor
            printf("Retrieving sequence #{ret_seq_count}, #{p_res.seq.id}, from NCBI at T+%02.0f:%02.0f:%02.0f\n",h,m,s)
            ncbi_blast_result = blaster.fetchTblastxResult(p_res)
            if(!ncbi_blast_result.nil?)
                ncbi_blast_results << ncbi_blast_result
            end
            put_results.delete(p_res)
            ret_seq_count += 1
        }
    end
}
if(put_results.length > 0)
    put_results.each {|rem_p_res|
        dt = Time.now - START_TIME
        h  = (dt / 3600).floor
        m  = ((dt % 3600) / 60).floor
        s  = ((dt % 3600) % 60).floor
        printf("Retrieving sequence #{ret_seq_count}, #{rem_p_res.seq.id}, from NCBI at T+%02.0f:%02.0f:%02.0f\n",h,m,s)
        ncbi_blast_result = blaster.fetchTblastxResult(rem_p_res)
        if(!ncbi_blast_result.nil?)
            ncbi_blast_results << ncbi_blast_result
        end
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
acc_num_groups = Hash.new
expr_sig_len = options[:ma_files].length
ncbi_blast_results.each {|ncbi_res|
    print "."
    $stdout.flush
    if(ncbi_res.hasAlignments?)
        acc_num = ncbi_res.bestAlignment.accession_num
        if(acc_num_groups[acc_num].nil?)
            acc_num_groups[acc_num] = AccessionNumGroup.new(acc_num,expr_sig_len,loghandl)
        end
        seq_id   = ncbi_res.sequence.id
        expr_sig = seq_hash[seq_id]
        acc_num_groups[acc_num].addRes(ncbi_res,expr_sig)
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
