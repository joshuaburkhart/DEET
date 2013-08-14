require_relative 'Sequence.rb'

class Alignment
    @sequence
    @accession_num
    @e_value
    def initialize(sequence,accession_num,e_value)
        if(!sequence.nil? && sequence.class == Sequence)
            @sequence = sequence
        else
            raise(ArgumentError,"ERROR: sequence must be of class 'Sequence'")
        end
        if(!accession_num.nil? && accession_num.class == String && accession_num.match(/^[0-9a-zA-Z._-]+$/))
            @accession_num = accession_num
        else
            raise(ArgumentError,"ERROR: invalid accession number")
        end
        if(!e_value.nil? && e_value.class == String && e_value.match(/^[0-9]+(.[0-9]+((e-)?[0-9]+)?)?$/))
            @e_value = e_value.to_f
        else
            raise(ArgumentError,"ERROR: invalid e value")
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
