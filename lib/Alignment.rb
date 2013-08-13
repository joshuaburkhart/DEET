class Alignment
    @sequence
    @accession_num
    @e_value
    def initialize(sequence,accession_num,e_value)
        @sequence = sequence
        @accession_num = accession_num
        @e_value = e_value
    end
end
