require_relative 'RgxLib'
class Sequence
    attr_accessor :id
    attr_accessor :bp_list

    attr_accessor :name
    attr_accessor :locus_tag
    attr_accessor :acc_num
    attr_accessor :expr_sig

    attr_accessor :ignored
    attr_accessor :orphan
    attr_accessor :paralog_num
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
        @name        = "N/A"
        @locus_tag   = "N/A"
        @expr_sig    = "N/A"
        @acc_num     = "N/A"
        @ignored     = true
        @orphan      = false
        @paralog_num = 0
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
    def getStatus
        status = nil
        if(@ignored)
            status = 0
        elsif(@orphan)
            status = 1
        elsif(@paralog_num == 0)
            status = 2
        elsif(@paralog_num > 0)
            status = 3
        end
        return status
    end
    def to_s
        row = "#{@id}~#{getStatus}~#{@name}~#{@locus_tag}~#{@acc_num}~#{@paralog_num}~#{@expr_sig}\n"
        return row
    end
end
