from Splice_Mod import *

class Display:
    def __init__(self,enzyme, rsq, min_l=0, max_l=1000, min_w=0, max_w=100000, mc=3): 
        self.enzyme = enzyme
        self.rsq = rsq
        self.min_l=min_l
        self.max_l=max_l
        self.min_w=min_w
        self.max_w=max_w
        self.mc = mc

        en = self.enzyme
        rs = self.rsq
        ml = self.min_l
        lm = self.max_l
        mw = self.min_w
        wm = self.max_w
        m_c = self.mc

        self.params = ['Enzyme:',en,'Legnth Range:',str(ml)+'-'+str(lm), 'Mass Range:', str(mw)+'-'+str(wm),'Missed_Cleavages:', str(m_c)]

        self.flag = 0
        if self.enzyme == 'Trypsin':
            Quer = Trypsin(rawseq=rs, min_l=ml, max_l=lm, min_w=mw, max_w=wm, mc=m_c)
            self.flag = 1           
        elif self.enzyme == 'Lys_C':
            Quer = Lys_C(rawseq=rs, min_l=ml, max_l=lm, min_w=mw, max_w=wm, mc=m_c)
            self.flag = 1
        elif self.enzyme == 'Lys_N':
            Quer = Lys_N(rawseq=rs, min_l=ml, max_l=lm, min_w=mw, max_w=wm, mc=m_c)
            self.flag = 1
        elif self.enzyme == 'CNBr':
            Quer = CNBr(rawseq=rs, min_l=ml, max_l=lm, min_w=mw, max_w=wm, mc=m_c)
            self.flag = 1
        elif self.enzyme == 'Arg_C':
            Quer = Arg_C(rawseq=rs, min_l=ml, max_l=lm, min_w=mw, max_w=wm, mc=m_c)
            self.flag = 1
        elif self.enzyme == 'Asp_N':
            Quer = Asp_N(rawseq=rs, min_l=ml, max_l=lm, min_w=mw, max_w=wm, mc=m_c)
            self.flag = 1
        self.table = Quer.LW_table()
