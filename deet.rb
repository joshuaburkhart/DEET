#!/usr/bin/ruby

require 'optparse'

require_relative 'lib/Sequence'
require_relative 'lib/Alignment'
require_relative 'lib/NCBIBlastResult'
require_relative 'lib/NCBIBlaster'
require_relative 'lib/AccessionNumGroup'
require_relative 'lib/FastaParser'
require_relative 'lib/MicroArrayHashBuilder'
require_relative 'lib/RgxLib'

EASTERN_OFFSET = 5
E_TIME = Time.now.localtime(EASTERN_OFFSET)
START_TIME = Time.now
run = false
if(!E_TIME.saturday? && !E_TIME.sunday? && !(E_TIME.hour > 21))
    puts "Due to the intensive load this program may put on NCBI servers, this program should only be run on weekends or between the hours of 9pm and 5am."
    puts "Choose One:"
    puts "w - wait for next eligible execution time, then run automatically"
    puts "a - abort execution and manually run this program later, or not"
    puts "r - run immediately, ignoring NCBI guidelines"
    answer = gets
    if(!answer.nil? && answer.class == String && answer.length == 1 && answer.match(/[war]/))
       if(answer == "w")
           num_minutes = 0
           t = Time.now.localtime(EASTERN_OFFSET)
           min_interval = 10
           puts "waiting..."
           while(!t.saturday? && !t.sunday? && !(t.hour > 21))
                 print "."
                 $stdout.flush
                 sleep(min_interval * 60)
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

MIN_SEQ_BP_LEN = 20

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

seqs = Array.new
options[:fasta_files].each {|fasta_file|
    puts "loading #{fasta_file}..."
    parser = FastaParser.new(fasta_file,MIN_SEQ_BP_LEN)
    parser.open
    while(next_seq = parser.nextSeq)
        seqs << next_seq
    end
    parser.close
}
puts "fasta files loaded."
puts "==================="

blaster = NCBIBlaster.new
put_results = Array.new
seqs.each {|seq|
    dt = Time.now - START_TIME
    h = (dt / 3600).floor
    m = ((dt % 3600) / 60).floor
    s = ((dt % 3600) % 60).floor
    printf("Time elapsed: %2.0f hours, %2.0f minutes, %2.0f seconds...\n",h,m,s)
    puts "submitting #{seq.id} to ncbi..."
    put_results << blaster.submitTblastxQuery(seq)
    sleep(3)
}
puts "sequences submitted to ncbi"
puts "==========================="

ncbi_blast_results = Array.new
put_results.each {|p_res|
    puts "retrieving #{p_res.seq.id} from ncbi..."
    ncbi_blast_results << blaster.fetchTblastxQuery(p_res)
    sleep(3)
}
puts "query results retreived"
puts "======================="

puts "hashing microarray data..."
seq_hash = MicroArrayHashBuilder.makeHash(options[:ma_files])
puts "microarray data hashed"
puts "======================"

puts "grouping sequences..."
acc_num_groups = Hash.new
expr_sig_len = options[:ma_files].length
ncbi_blast_results.each {|ncbi_res|
    print "."
    $stdout.flush
    acc_num = ncbi_res.bestAlignment.accession_num
    if(acc_num_groups[acc_num].nil?)
        acc_num_groups[acc_num] = AccessionNumGroup.new(acc_num,expr_sig_len)
    end
    seq_id = ncbi_res.sequence.id
    expr_sig = seq_hash[seq_id]
    acc_num_groups[acc_num].addRes(acc_num,expr_sig)
}
puts
puts "sequences grouped by accession number and expression signature"
puts "=============================================================="

puts "writing results to hard disk..."
outhandl = File.open("deet.#{Time.now.to_i}.result","w")
acc_num_groups.values.each {|acc_num_group|
    print "."
    $stdout.flush
    outhandl.puts "-----------------"
    outhandl.print(acc_num_group.to_s)
    outhandl.puts "-----------------"
}
puts
outhandl.close
puts "done."
