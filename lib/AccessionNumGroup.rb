require_relative 'RgxLib'
require_relative 'Sequence'

class AccessionNumGroup
    @acc_num
    @expr_sig_len
    @paralogs = Hash.new
    def initialize(acc_num,expr_sig_len)
        if(!acc_num.nil? && acc_num.class == String && acc_num.match(RgxLib::ALGN_ACC_NUM))
            @acc_num = acc_num
            if(!expr_sig_len.nil? && expr_sig_len.class == Fixnum && expr_sig_len > 0)
                @expr_sig_len = expr_sig_len
            else
                raise(ArgumentError,"ERROR: expr_sig_len '#{expr_sig_len}' not a valid length")
            end
        else
            raise(ArgumentError,"ERROR: acc_num '#{acc_num}' not a valid accession number")
        end
    end
    def addRes(expr_sig,res)
        if(!expr_sig.nil? && expr_sig.class == String && expr_sig.match(RgxLib::ACCG_EXPR_SIG))
            if(!res.nil? && res.class == NCBIBlastResult)
                if(@paralogs[expr_sig].nil?)
                    @paralogs[expr_sig] = Array.new
                end
                @paralogs[expr_sig] << res
            else
                raise(ArgumentError,"ERROR: res not valid")
            end
        else
            raise(ArgumentError,"ERROR: expr_sig '#{expr_sig}' not a valid expression signature")
        end
    end
    def getRepresentativeSeq(expr_sig)
        if(!expr_sig.nil? && expr_sig.calss = String \
           && expr_sig.length == @expr_sig_len \
           && expr_sig.match(RgxLib::ACCG_EXPR_SIG))
            return @paralogs[expr_sig]
        else
            raise(ArgumentError,"ERROR: expr_sig '#{expr_sig}' not a valid expression signature")
        end
        def getParalogExprSigs
            gene_sig = getGeneExprSig
            return @paralogs.to_a - [[gene_sig,@paralogs[gene_sig]]]
        end
        def getGeneExprSig
            lowest_e_val = nil
            gene_sig = nil
            @paralogs.to_a.each {|paralog|
                paralog.each {|ncbi_result|
                    cur_e_val = ncbi_result.bestAlignment.e_value
                    if(gene.nil? || cur_e_val < lowest_e_val)
                        lowest_e_val = cur_e_val
                        gene_sig = @paralogs.key(paralog)
                    end
                }
            }
            return gene_sig
        end
    end
end
