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

EASTERN_OFFSET = (-5 * 3600)
E_TIME = Time.now.localtime(EASTERN_OFFSET)
START_TIME = Time.now
MIN_SEQ_BP_LEN = 100
OUT_DIR_NAME = "output/"
%x(mkdir -p #{OUT_DIR_NAME})
result_filename = "#{OUT_DIR_NAME}#{START_TIME.to_i}.result"
log_filename = "#{OUT_DIR_NAME}#{START_TIME.to_i}.log"

if(!E_TIME.saturday? && !E_TIME.sunday? && !(E_TIME.hour > 21))
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
            t = Time.now.localtime(EASTERN_OFFSET)
            puts "waiting..."
            while(!t.saturday? && !t.sunday? && !(t.hour > 21))
                print "."
                $stdout.flush
                sleep(600) #10 minutes
                t = Time.now.localtime(EASTERN_OFFSET)
            end         
            puts "beginning execution after #{num_minutes} minute wait..."
            run = true
        elsif(answer == "a")
            puts "aborted."
            exit(0)
        elsif(answer == "r")
            puts "good luck..."
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
        puts "accepted -f '#{options[:fasta_files]}' as argument."
    }
    options[:ma_files] = nil
    opts.on('-m','--mas MA_FILE1,MA_FILE2',Array,'MA files') { |ma_files|
        options[:ma_files] = ma_files
        puts "accepted -m '#{options[:ma_files]}' as argument."
    }
}
optparse.parse!

if(options[:fasta_files].nil?)
    raise(ArgumentError,"ERROR: no FASTA files supplied")
elsif(options[:ma_files].nil?)
    raise(ArgumentError,"ERROR: no MA files supplied")
end

loghandl = File.open(log_filename,"w")
puts "hashing microarray data..."
seq_hash = MicroArrayHashBuilder.makeHash(*options[:ma_files],loghandl)
puts "microarray data hashed"
puts "======================"

seqs = Set.new
options[:fasta_files].each {|fasta_file|
    puts "loading #{fasta_file}..."
    parser = FastaParser.new(fasta_file,MIN_SEQ_BP_LEN,loghandl)
    parser.open
    while(next_seq = parser.nextSeq)
        seqs << next_seq
        print "."
        $stdout.flush
    end
    parser.close
}
puts
puts "fasta files loaded. (#{seqs.length} total sequences)"
puts "=================================================="

puts "querying ncbi..."
blaster = NCBIBlaster.new(loghandl)
put_results = Set.new
ncbi_blast_results = Set.new
ret_seq_count = 0
seqs.each_with_index {|seq,i|
    dt = Time.now - START_TIME
    h = (dt / 3600).floor
    m = ((dt % 3600) / 60).floor
    s = ((dt % 3600) % 60).floor
    printf("submitting sequence #{i}, #{seq.id}, to ncbi at T+%02.0f:%02.0f:%02.0f\n",h,m,s)
    put_results << blaster.submitTblastxQuery(seq)
    if(i % 100 == 99)
        put_results.each {|p_res|
            dt = Time.now - START_TIME
            h = (dt / 3600).floor
            m = ((dt % 3600) / 60).floor
            s = ((dt % 3600) % 60).floor
            printf("retrieving sequence #{ret_seq_count}, #{p_res.seq.id}, from ncbi at T+%02.0f:%02.0f:%02.0f\n",h,m,s)
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
        h = (dt / 3600).floor
        m = ((dt % 3600) / 60).floor
        s = ((dt % 3600) % 60).floor
        printf("retrieving sequence #{ret_seq_count}, #{rem_p_res.seq.id}, from ncbi at T+%02.0f:%02.0f:%02.0f\n",h,m,s)
        ncbi_blast_result = blaster.fetchTblastxResult(rem_p_res)
        if(!ncbi_blast_result.nil?)
            ncbi_blast_results << ncbi_blast_result
        end
        put_results.delete(rem_p_res)
        ret_seq_count += 1
    }
end
puts "query results retreived"
puts "======================="

puts "grouping sequences..."
acc_num_groups = Hash.new
expr_sig_len = options[:ma_files].length
ncbi_blast_results.each {|ncbi_res|
    print "."
    $stdout.flush
    if(ncbi_res.hasAlignments?)
        acc_num = ncbi_res.bestAlignment.accession_num
        if(acc_num_groups[acc_num].nil?)
            acc_num_groups[acc_num] = AccessionNumGroup.new(acc_num,expr_sig_len)
        end
        seq_id = ncbi_res.sequence.id
        expr_sig = seq_hash[seq_id]
        acc_num_groups[acc_num].addRes(ncbi_res,expr_sig)
    end
}
puts
puts "sequences grouped by accession number and expression signature"
puts "=============================================================="

puts "writing results to #{result_filename}..."
resulthandl = File.open(result_filename,"w")
acc_num_groups.values.each {|acc_num_group|
    print "."
    $stdout.flush
    resulthandl.puts(acc_num_group.to_s)
    resulthandl.puts
}
puts
resulthandl.close
loghandl.close
puts "done."
