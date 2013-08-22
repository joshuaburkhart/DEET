module RgxLib
    SEQ_ID = /^([^\s\\n]+)$/
    SEQ_BP_LIST = /^([atcgATCG]+)$/    
    ALGN_ACC_NUM = /^[0-9a-zA-Z._-]+$/
    ALGN_E_VAL = /^([0-9.]+(e-)?([0-9]+)?)$/
    BLST_NO_MATCH = /No significant similarity found/
    BLST_ACCN_GRAB = /\|([^\|]+)\|/
    BLST_E_VAL_GRAB = /.+?[a-z]+\|[0-9a-zA-Z._-]+\|\s+[\S ]+\s+[0-9.]+\s+([0-9.]+(e-)?([0-9]+)?)\s+[0-9]+$/m
    BLST_RID_GRAB = /RID = ([0-9A-Z-]+)/
    BLST_RTOE_GRAB = /RTOE = ([0-9]+)/
    BLST_WAIT = /Status=WAITING/
    BLST_READY = /Status=READY/
    BLST_HEADER_GRAB = /(<p>.*?Query=)/m
    ACCG_EXPR_SIG = /[X10]+/
    MAHB_HEADER_CHECK = /^gene\s+M\s+A\s+t\s+pvalue\s+qvalue\s+B$/
    MAHB_SEQ_ID_GRAB = /^([FC][0-9A-Z]+).+$/
    MAHB_Q_VAL_GRAB = /^[FC][0-9A-Z]+\s+-?[0-9.]+(e-)?([0-9]+)?\s+-?[0-9.]+(e-)?([0-9]+)?\s+-?[0-9.]+(e-)?([0-9]+)?\s+-?[0-9.]+(e-)?([0-9]+)?\s+(-?[0-9.]+(e-)?([0-9]+)?).+$/
    MAHB_M_VAL_GRAB = /^[FC][0-9A-Z]+\s+(-?[0-9.]+(e-)?([0-9]+)?).+$/
end
