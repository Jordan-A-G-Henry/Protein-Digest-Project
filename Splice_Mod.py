#The purpose of this module is to accept a protein sequence, and output a table of information.

class Splice_Prot:
    def __init__(self,rawseq, min_l=0, max_l=1000, min_w=0, max_w=100000, mc=0): 
        self.rawseq = rawseq
        self.min_l=min_l
        self.max_l=max_l
        self.min_w=min_w
        self.max_w=max_w
        #A placeholder. Each child class builds up its own list of splice sites.
        self.splice_list = [0]

    #Generates a list of start and end points for each splice.
    #Accounts for the missed cleavages condition to create all possible splices based on this data.
    #Together with the sequence itself, this is all the information needed to generate the table.
    def full_list(self):
        #Add end of sequence to splice list.
        self.splice_list.append(len(self.rawseq))            
        tuple_list = []
        for i in range(0,len(self.splice_list)-1):
            for j in range (0,self.mc+1):
                if not ((i+1+j)+1) > len(self.splice_list):
                    #Build up the entry (as a list).
                    entry = []                                 
                    entry.append(self.splice_list[i])
                    entry.append(self.splice_list[i+1+j])
                    entry.append(j)                         
                    #Add it to the list (of lists).
                    tuple_list.append(entry)           
        return tuple_list

    #Create the actual table.
    def LW_table(self):
        splice_table = [['Start-Stop','Legnth'
                       ,'Mol_Weight','Missed_Cleavages','Left-Of-Seq'
                       ,'Peptide','Right-Of-Seq']]

        full_list = self.full_list()
        #For some reason, duplicates are always added to the list.
        clean_list = []
        for i in full_list:
            if not i in clean_list:
                clean_list.append(i)
        full_list = clean_list

        for i in range(0,(len(full_list)-1)):
            entry = []
            #1.Start-stop position.
            entry.append((str(full_list[i][0]+1))+'-'+(str(full_list[i][1])))
            #Bonus: Legnth
            entry.append(full_list[i][1] - full_list[i][0])
            #2.Molecular weight.
            mw_tot = 19.02
            mw = {'A': 71.04, 'C': 103.01, 'D': 115.03, 'E': 129.04, 'F': 147.07,
          	  'G': 57.02, 'H': 137.06, 'I': 113.08, 'K': 128.09, 'L': 113.08,
          	  'M': 131.04, 'N': 114.04, 'P': 97.05, 'Q': 128.06, 'R': 156.10,
          	  'S': 87.03, 'T': 101.05, 'V': 99.07, 'W': 186.08, 'Y': 163.06 }
            for n in range(full_list[i][0],full_list[i][1]):
               mw_tot+= mw[self.rawseq[n]]
               mw_tot = round(mw_tot,2)                 
            entry.append(mw_tot)
            #3.Missed cleavages. Taken directly from previous list.
            entry.append(full_list[i][2])          
            #4. Left-of-Sequence.
            l = (full_list[i][0])-1
            if l < 0:
                entry.append('_')
            else:
                entry.append(self.rawseq[l])
            #5. Actual-Sequence.
            entry.append(self.rawseq[full_list[i][0]:full_list[i][1]])
            #6. Right-of-Squence.
            r = full_list[i][1]+1
            try:
                entry.append(self.rawseq[r])
            except IndexError:
                entry.append('*')

            
            #Finally: For each entry, use length and weight conditions to decide whether to add it to the table.
            if entry[1]>self.min_l and entry[1]<self.max_l:
                if entry[2]>self.min_w and entry[2]<self.max_w:
                    splice_table.append(entry)


        return (splice_table)



#All the child classes now inherit the above methods. Because they each feed into them a different list of splice sites,
#each determined according to the logic of the actual enzyme, the result is a different final table.
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




