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

MIN_SEQ_BP_LEN = 20

options = {}
optparse = OptionParser.new { |opts|
    opts.banner = <<-EOS
Usage: ruby deet.rb -f <fasta file 1> ... <fasta file n> -m <ma file 1> ... <ma file n>

Example:
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
    opts.on('-m','--mas MA_FILE1,MA_FILE2',Array,'MA files') { |ma_files|
        options[:ma_files] = ma_files
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
    puts "submitting #{seq.id} to ncbi..."
    put_results << blaster.submitTblastxQuery(seq)
}
puts "sequences submitted to ncbi"
puts "==========================="

ncbi_blast_results = Array.new
put_results.each {|p_res|
    puts "retrieving #{p_res.seq.id} from ncbi..."
    ncbi_blast_results << blaster.fetchTblastxQuery(p_res)
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
