class NCBIBlastResult
    @sequence
    @alignments
    def initialize(sequence)
        @sequence = sequence
        @alignments = Array.new
    end
    def addAlignment(alignment)
        @alignments << alignment
    end
    def alignmentCount()
        return @alignments.length
    end
    def bestAlignment()
        @alignments.sort {|i,j| i.e_value <=> j.e_value}[0]
    end
end
