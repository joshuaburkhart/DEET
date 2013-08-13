class Sequence
    @id
    @bp_list
    def initialize(id,bp_list)
        @id = id
        @bp_list = bp_list
    end
    def nil?()
        return @id.nil? || @bp_list.nil?
    end
end
