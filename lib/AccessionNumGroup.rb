require_relative 'RgxLib'
require_relative 'Sequence'

class AccessionNumGroup
    @acc_num
    @expr_sig_len
    @paralogs
    def initialize(acc_num,expr_sig_len)
        if(!acc_num.nil? && acc_num.class == String && acc_num.match(RgxLib::ALGN_ACC_NUM))
            if(!expr_sig_len.nil? && expr_sig_len.class == Fixnum && expr_sig_len > 0)
                @acc_num = acc_num
                @expr_sig_len = expr_sig_len
                @paralogs = Hash.new
            else
                raise(ArgumentError,"ERROR: expr_sig_len '#{expr_sig_len}' not a valid length")
            end
        else
            raise(ArgumentError,"ERROR: acc_num '#{acc_num}' not a valid accession number")
        end
    end
    def addRes(res,expr_sig)
        if(!expr_sig.nil? && expr_sig.class == String && expr_sig.match(RgxLib::ACCG_EXPR_SIG))
            if(!res.nil? && res.class == NCBIBlastResult && \
               !res.bestAlignment.nil? && res.bestAlignment.accession_num == @acc_num)
                if(@paralogs[expr_sig].nil?)
                    @paralogs[expr_sig] = Array.new
                end
                @paralogs[expr_sig] << res
            else
                raise(ArgumentError,"ERROR: res '#{res}' not valid")
            end
        else
            raise(ArgumentError,"ERROR: expr_sig '#{expr_sig}' not a valid expression signature")
        end
    end
    def getRepresentativeSeq(expr_sig)
        if(!expr_sig.nil? && expr_sig.class == String \
           && expr_sig.length == @expr_sig_len \
           && expr_sig.match(RgxLib::ACCG_EXPR_SIG))
            return getRepresentativeRes(expr_sig).sequence
        else
            raise(ArgumentError,"ERROR: expr_sig '#{expr_sig}' not a valid expression signature")
        end
    end
    def getRepresentativeRes(expr_sig)
        if(!expr_sig.nil? && expr_sig.class == String \
           && expr_sig.length == @expr_sig_len \
           && expr_sig.match(RgxLib::ACCG_EXPR_SIG))
            return longestRes(@paralogs[expr_sig])
        else
            raise(ArgumentError,"ERROR: expr_sig '#{expr_sig}' not a valid expression signature")
        end
    end
    def getParalogExprSigs
        gene_sig = getGeneExprSig
        return @paralogs.keys - [gene_sig]
    end
    def getGeneExprSig
        lowest_e_val = Float::INFINITY
        gene_sig = nil
        @paralogs.values.each {|paralog|
            paralog.each {|ncbi_result|
                cur_e_val = ncbi_result.hasAlignments? ? ncbi_result.bestAlignment.e_value : Float::INFINITY
                if(gene_sig.nil? || cur_e_val < lowest_e_val)
                    lowest_e_val = cur_e_val
                    gene_sig = @paralogs.key(paralog)
                end
            }
        }
        return gene_sig
    end
    def longestRes(ncbi_res_ary)
        if(!ncbi_res_ary.nil? && ncbi_res_ary.class == Array && ncbi_res_ary.length > 0)
            res = ncbi_res_ary.sort {|i,j| j.sequence.bp_list.length <=> i.sequence.bp_list.length }[0]
            return res
        else
            raise(ArgumentError,"ERROR: ncbi_res_ary '#{ncbi_res_ary}' not a valid array")
        end
    end
    def to_s
        gene_expr_sig = getGeneExprSig
        best_gene_rep = getRepresentativeSeq(gene_expr_sig)
        best_gene_ncbi_res = getRepresentativeRes(gene_expr_sig)
        string_rep =<<EOS
GENE: #{@acc_num}
======#{'='*@acc_num.length}
\tSIGNATURE: #{gene_expr_sig}
\tBEST REPRESENTATIVE SEQUENCE ID: #{best_gene_rep.id}
\tBEST REPRESENTATIVE SEQUENCE LENGTH: #{best_gene_rep.bp_list.length}
\tBEST REPRESENTATIVE SEQUENCE E VALUE: #{best_gene_ncbi_res.bestAlignment.e_value}
EOS
        other_gene_res = @paralogs[gene_expr_sig] - [best_gene_ncbi_res]
        other_gene_res.each_with_index {|res,r|
            other_gene_rep =<<EOS
\t\tREPRESENTATIVE #{r + 2} SEQUENCE ID: #{res.sequence.id}
\t\tREPRESENTATIVE #{r + 2} SEQUENCE LENGTH: #{res.sequence.bp_list.length}
\t\tREPRESENTATIVE #{r + 2} SEQUENCE E VALUE: #{res.bestAlignment.e_value}
EOS
            string_rep += other_gene_rep
        }
        getParalogExprSigs.each_with_index {|psig,i|
            paralog_expr_sig = psig
            best_paralog_rep = getRepresentativeSeq(psig)
            best_paralog_ncbi_res = getRepresentativeRes(psig)
            p_rep =<<EOS
\t#{@acc_num} PARALOG #{i + 1}
\t#{'-'*@acc_num.length}----------
\t\tSIGNATURE: #{psig}
\t\tBEST REPRESENTATIVE SEQUENCE ID: #{best_paralog_rep.id}
\t\tBEST REPRESENTATIVE SEQUENCE LENGTH: #{best_paralog_rep.bp_list.length}
\t\tBEST REPRESENTATIVE SEQUENCE E VALUE: #{best_paralog_ncbi_res.bestAlignment.e_value}
EOS
            other_paralog_res = @paralogs[psig] - [best_paralog_ncbi_res]
            other_paralog_res.each_with_index {|pres,j|
                other_paralog_rep =<<EOS
\t\t\tREPRESENTATIVE #{j + 2} SEQUENCE ID: #{pres.sequence.id} 
\t\t\tREPRESENTATIVE #{j + 2} SEQUENCE LENGTH: #{pres.sequence.bp_list.length} 
\t\t\tREPRESENTATIVE #{j + 2} SEQUENCE E VALUE: #{pres.bestAlignment.e_value} 
EOS
                p_rep += other_paralog_rep
            }
            string_rep += p_rep
        }
        return string_rep
    end
end
