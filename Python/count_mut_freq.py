#!/usr/local/bin/ python2.7
# encoding: utf-8

import os,sys

# HCV_RAMS FOR BOC, TVR
NS3_RAMS = {'V36A':'','V36I':'','V36M':'','V36G':'','V36C':'','V36L':'','Q41R':'','F43C':'','F43S':'','T54A':'','T54S':'','T54G':'','T54C':'','V55A':'','V55I':'','V107I':'','I132V':'',
            'R155G':'','R155I':'','R155K':'','R155M':'','R155Q':'','R155T':'','R155P':'','R155LL':'','R155S':'','A156S':'','A156T':'','A156V':'','A156G':'','A156F':'','A156N':'',
            'A156I':'','V158I':'','D168N':'','D168Y':'','I170A':'','I170F':'','I170T':'','I170L':'','M175L':''}

ram_count = {}
co = {}

def checkSingleMut(mut):
    if mut != 'None' and NS3_RAMS.has_key(mut):
        if ram_count.has_key(mut):
            ram_count[mut] = ram_count[mut] + 1
        else:
            ram_count[mut] = 1
        return 'yes'
    return 'no'
    
    
def countMuts(mut_sum):
    if str(mut_sum) !='None':
        ramstring = ''
        muts = mut_sum.replace('"','').split(', ')
        for m in muts:
            # Check for single mutants first
            if m.find('/') == -1: 
                ram = checkSingleMut(m)
                if ram =='yes': ramstring = ramstring+m+','
            
            # Have multiple mutants
            else:
                muts = m.split('/')   # ['I50A','C','H']
                mut_count = len(muts)
                
                # check first mutation in list
                if muts[0][0] != muts[0][-1]: 
                    ram = checkSingleMut(muts[0])
                    if ram =='yes': 
                        ramstring = ramstring + muts[0]+','
                
                # Generate the rest of mutations in list and check them.
                for i in range(1,len(muts)):
                    m2 = muts[0][:-1] + muts[i] # automatically create I50C, I50H etc
                    ram = checkSingleMut(m2)
                    if ram == 'yes': ramstring = ramstring+m2+','
                

        if ramstring != '':
            ramstring = ramstring[:-1]
            if co.has_key(ramstring): co[ramstring] = co[ramstring] + 1
            else: co[ramstring] = 1

    return


def main():
    MUTFILE = open('500_samples.txt','r')
    OUTFILE = open('500_samples_counted.txt','w')   
    
    for line in MUTFILE:
        vals = line.split('\t')
        ns3_muts = vals[5]

        countMuts(ns3_muts)
    
    rk = ram_count.keys()
    rk.sort()
    for r in rk:
        OUTFILE.write(r+'\t'+str(ram_count[r])+'\n')
    OUTFILE.write("\n\nCo-occurance\n")
    ck = co.keys()
    ck.sort()
    for k in ck:
        OUTFILE.write(k+'\t'+str(co[k])+'\n')
    OUTFILE.close()
    MUTFILE.close()


if __name__ == '__main__':
   main()