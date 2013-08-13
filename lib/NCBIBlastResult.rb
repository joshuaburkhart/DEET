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
        return @alignments[0].e_value > @alignments[1].e_value ? @alignments[1] : @alignments[0]
    end
end
