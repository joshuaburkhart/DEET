require_relative 'RgxLib'
class Sequence
    attr_accessor :id
    attr_accessor :bp_list
    def initialize(id,bp_list)
        if(id.match(RgxLib::SEQ_ID))
           @id = $1
        else
            raise(ArgumentError, "Sequence could not be initialized with non-compliant id '#{id}'.")
        end
        if(bp_list.match(RgxLib::SEQ_BP_LIST))
        @bp_list = $1
        else
            raise(ArgumentError, "Sequence could not be initialized with non-compliant bp_list '#{bp_list}'")
        end
    end
    def nil?
        return @id.nil? || @bp_list.nil?
    end
    def eql?(other)
        return (!@id.nil? && other.class == Sequence && @id == other.id)
    end
    def hash
        return (@id.hash)
    end
end
