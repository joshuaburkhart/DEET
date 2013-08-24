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
    def addRes(expr_sig,res)
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
            return longestSeq(@paralogs[expr_sig])
        else
            raise(ArgumentError,"ERROR: expr_sig '#{expr_sig}' not a valid expression signature")
        end
    end
    def getParalogExprSigs
        gene_sig = getGeneExprSig
        return @paralogs.values - [@paralogs[gene_sig]]
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
    def longestSeq(ncbi_res_ary)
        if(!ncbi_res_ary.nil? && ncbi_res_ary.class == Array && ncbi_res_ary.length > 0)
            res = ncbi_res_ary.sort {|i,j| j.sequence.bp_list.length <=> i.sequence.bp_list.length }[0]
            return res.sequence
        else
            raise(ArgumentError,"ERROR: ncbi_res_ary '#{ncbi_res_ary}' not a valid array")
        end
    end
    def to_s
        string_rep =<<EOS
GENE ACCESSION NUMBER:     #{@acc_num}
GENE EXPRESSION SIGNATURE: #{getGeneExprSig}
GENE REPRESENTATIVE:       #{getRepresentativeSeq(getGeneExprSig)}
EOS
        getParalogExprSigs.each_with_index {|psig,i|
        p_rep =<<EOS
\tPARALOG_#{i} EXPRESSION SIGNATURE: #{psig}
\tPARALOG_#{i} REPRESENTATIVE:       #{getRepresentativeSeq(psig)}
EOS
        string_rep += p_rep
        }
        return string_rep
    end
end
