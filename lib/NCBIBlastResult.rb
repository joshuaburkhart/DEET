require_relative 'Sequence.rb'

class NCBIBlastResult
    attr_accessor :valid
    @sequence
    @alignments
    def initialize(sequence)
        if(!sequence.nil? && sequence.class == Sequence)
            @sequence = sequence
            @valid = true
        else
            raise(ArgumentError,"ERROR: sequence must be an instance of 'Sequence' not '#{sequence.class}'")
        end
        @alignments = Array.new
    end
    def alignments
        return @alignments
    end
    def sequence
        return @sequence
    end
    def addAlignment(alignment)
        @alignments << alignment
    end
    def hasAlignments?
        return (alignmentCount > 0)
    end
    def alignmentCount
        return @alignments.length
    end
    def bestAlignment
        @alignments.sort {|i,j| i.e_value <=> j.e_value}[0]
    end
    def eql?(other)
        #NCBIResults can be equal to a sequence
        return (!@sequence.nil? && !@sequence.id.nil? && other.class == Sequence && @id == other.id)
    end 
end
