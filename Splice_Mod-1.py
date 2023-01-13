#The purpose of this module is to accept a protein sequence, and output a table of information.



#[1]

class Splice_Prot:
    def __init__(self,rawseq, min_l=0, max_l=1000, min_w=0, max_w=100000, mc=0): 
        #Initialize user inputs. Not the default values.
        self.rawseq = rawseq
        self.min_l=min_l
        self.max_l=max_l
        self.min_w=min_w
        self.max_w=max_w
        self.mc = mc
#[2]
        self.splice_list = [0,10,20,30,40]

    def table(self):
        splice_table = [['Start-Stop'
                       ,'Mol_Weight','Missed_Cleavages','Left-Of-Seq'
                       ,'Peptide','Right-Of-Seq']]
#[3]
        self.splice_list.append(len(self.rawseq))
        for i in range(0,(len(self.splice_list)-1)):     
            j = i+1
#[4]
            entry = []
#[4.1]
            entry.append((str(self.splice_list[i]+1))+'-'+(str(self.splice_list[j]))) 
#[4.2]   
            mw_tot = 0.
            mw = {'A': 71.04, 'C': 103.01, 'D': 115.03, 'E': 129.04, 'F': 147.07,
          	  'G': 57.02, 'H': 137.06, 'I': 113.08, 'K': 128.09, 'L': 113.08,
          	  'M': 131.04, 'N': 114.04, 'P': 97.05, 'Q': 128.06, 'R': 156.10,
          	  'S': 87.03, 'T': 101.05, 'V': 99.07, 'W': 186.08, 'Y': 163.06 }
            for n in range(self.splice_list[i],self.splice_list[j]):
                try:
                    mw_tot+= mw[self.rawseq[n]]
                except KeyError:
                    pass      
            entry.append(mw_tot)
#[4.3]
            entry.append(0)
#[4.4]
            l = self.splice_list[i]-1
#[5]
            if i == 0:
                entry.append('_')
            else:
                entry.append(self.rawseq[l])
#[4.5]
            entry.append(self.rawseq[self.splice_list[i]:self.splice_list[j]])
#[4.6]
            r = self.splice_list[j]+1
            try:
                entry.append(self.rawseq[r])
            except IndexError:
                entry.append('*')            
#[6]
            splice_table.append(entry)
        return (splice_table)


    
#[7]
    def MC_table(self):
        In_table = self.table()
#[8]
        splice_table = [['Start-Stop'
                       ,'Mol_Weight','Missed_Cleavages','Left-Of-Seq'
                       ,'Peptide','Right-Of-Seq']]
        for i in range (1,len(In_table)):
            type(In_table[i][1])
            if In_table[i][1] > 0:
                splice_table.append(In_table[i])

#[9]
        MC_table = [['Start-Stop'
                       ,'Mol_Weight','Missed_Cleavages','Left-Of-Seq'
                       ,'Peptide','Right-Of-Seq']]
        MC_table.append(splice_table[1])
#[10]
        for i in range (2,len(splice_table)):  
            MC_table.append(splice_table[i])

#[11]
            add_entry = []
#[11.1]  
            add_entry.append('('+splice_table[i-1][0]+'),('+splice_table[i][0]+')')
#[11.2]       
            add_entry.append(splice_table[i-1][1]+splice_table[i][1])
#[11.3]       
            add_entry.append(1)
#[11.4]        
            add_entry.append(splice_table[i-1][3])
#[11.5]          
            add_entry.append(splice_table[i-1][4]+splice_table[i][4])
#[11.6]
            add_entry.append(splice_table[i][5])
#[11.7]            
            if self.mc > 0:
                MC_table.append(add_entry)                       


#[12]
            add_matrix = []
            add_matrix.append(add_entry)
            for j in range(1, self.mc):
                if i > (j+1):                           
                    add_matrix.append([0,0,0,0,0,0])
                    add_matrix[j][0] = '('+splice_table[i-j-1][0]+'),' + add_matrix[j-1][0]
                    add_matrix[j][1] = add_matrix[j-1][1] + splice_table[i-j-1][1]
                    add_matrix[j][2] = add_matrix[j-1][2] + 1
                    add_matrix[j][3] = splice_table[i-j-1][3]
                    add_matrix[j][4] = splice_table[i-j-1][4] + add_matrix[j-1][4]
                    add_matrix[j][5] = add_matrix[j-1][5] 

                    MC_table.append(add_matrix[j])
            
        return MC_table

#[13]
    def LW_table(self):
        In_table = self.MC_table()

        splice_table = [['Start-Stop'
                       ,'Mol_Weight','Missed_Cleavages','Left-Of-Seq'
                       ,'Peptide','Right-Of-Seq']]    

        for i in range (1,len(In_table)):
            type(In_table[i][1])
            if In_table[i][1] > 0:
                splice_table.append(In_table[i])

        LW_table = [['Start-Stop'
                       ,'Mol_Weight','Missed_Cleavages','Left-Of-Seq'
                       ,'Peptide','Right-Of-Seq']]

        ml = self.min_l
        lm = self.max_l
        mw = self.min_w
        wm = self.max_w
        for i in range (1,len(In_table)):
            if len(In_table[i][4]) > ml and len(In_table[i][4])< lm:
                if In_table[i][1] > mw and In_table[i][1] < wm:
                    LW_table.append(In_table[i])

        return LW_table


#[14]
class Trypsin(Splice_Prot):
    def __init__(self,rawseq, min_l=0, max_l=1000, min_w=0, max_w=100000, mc=0): 
        self.rawseq = rawseq
        self.min_l=min_l
        self.max_l=max_l
        self.min_w=min_w
        self.max_w=max_w
        self.mc = mc
        self.splice_list = [0]         
        for n in range(0,len(self.rawseq)):
            if self.rawseq[n] == 'K' or self.rawseq[n] == 'R':
                if self.rawseq[n+1] != 'P':
                    self.splice_list.append(int(n+1))

class Lys_C(Splice_Prot):
    def __init__(self,rawseq, min_l=0, max_l=1000, min_w=0, max_w=100000, mc=0): 
        self.rawseq = rawseq
        self.min_l=min_l
        self.max_l=max_l
        self.min_w=min_w
        self.max_w=max_w
        self.mc = mc
        self.splice_list = [0]         
        for n in range(0,len(self.rawseq)):
            if self.rawseq[n] == 'K':
                    self.splice_list.append(int(n+1))

class Lys_N(Splice_Prot):
    def __init__(self,rawseq, min_l=0, max_l=1000, min_w=0, max_w=100000, mc=0): 
        self.rawseq = rawseq
        self.min_l=min_l
        self.max_l=max_l
        self.min_w=min_w
        self.max_w=max_w
        self.mc = mc
        self.splice_list = [0]         
        for n in range(0,len(self.rawseq)):
            if self.rawseq[n] == 'K':
                    self.splice_list.append(int(n))

class CNBr(Splice_Prot):
    def __init__(self,rawseq, min_l=0, max_l=1000, min_w=0, max_w=100000, mc=0): 
        self.rawseq = rawseq
        self.min_l=min_l
        self.max_l=max_l
        self.min_w=min_w
        self.max_w=max_w
        self.mc = mc
        self.splice_list = [0]         
        for n in range(0,len(self.rawseq)):
            if self.rawseq[n] == 'M':
                    self.splice_list.append(int(n+1))

class Arg_C(Splice_Prot):
    def __init__(self,rawseq, min_l=0, max_l=1000, min_w=0, max_w=100000, mc=0): 
        self.rawseq = rawseq
        self.min_l=min_l
        self.max_l=max_l
        self.min_w=min_w
        self.max_w=max_w
        self.mc = mc
        self.splice_list = [0]         
        for n in range(0,len(self.rawseq)):
            if self.rawseq[n] == 'R':
                if self.rawseq[n+1] == 'R':
                    self.splice_list.append(int(n+1))

class Asp_N(Splice_Prot):
    def __init__(self,rawseq, min_l=0, max_l=1000, min_w=0, max_w=100000, mc=0): 
        self.rawseq = rawseq
        self.min_l=min_l
        self.max_l=max_l
        self.min_w=min_w
        self.max_w=max_w
        self.mc = mc
        self.splice_list = [0]         
        for n in range(0,len(self.rawseq)):
            if self.rawseq[n] == 'D':
                    self.splice_list.append(int(n))
    


