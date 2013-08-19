require_relative 'RgxLib'
require_relative 'Sequence.rb'

class Alignment
    @sequence
    @accession_num
    @e_value
    def initialize(sequence,accession_num,e_value)
        if(!sequence.nil? && sequence.class == Sequence)
            @sequence = sequence
        else
            raise(ArgumentError,"ERROR: sequence must be of class 'Sequence' but was '#{sequence.class}'")
        end
        if(!accession_num.nil? && accession_num.class == String && accession_num.match(RgxLib::ALGN_ACC_NUM))
            @accession_num = accession_num
        else
            raise(ArgumentError,"ERROR: invalid accession number '#{accession_num}'")
        end
        if(!e_value.nil? && e_value.class == String && e_value.match(RgxLib::ALGN_E_VAL))
            @e_value = e_value.to_f
        else
            raise(ArgumentError,"ERROR: invalid e value '#{e_value}'")
        end
    end
    def sequence
        return @sequence
    end
    def accession_num
        return @accession_num
    end
    def e_value
        return @e_value
    end
end
