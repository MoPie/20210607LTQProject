#!/usr/bin/env python

import sys, optparse , time , numpy , re, resource

start_time = time.time()
##############################################
#    DMR calling using Fisher exact test     #
#       by Wanlu Liu                         #
#        from Jacobsen Lab, UCLA             #
##############################################
#Algorithm for DMRcalling in this script.
############################################################################################################################################
#INPUT                                                                                                                                     #
#1. Take input from BSMAP(from methratio.py) output (xxx_methratio.txt)                                                                    #
#chr     pos     strand  context ratio   eff_CT_count    C_count CT_count        rev_G_count     rev_GA_count    CI_lower        CI_upper  #
#chr1    15      +       AACCC   0.000   1.00    0       1       8       8       0.000   0.793                                             #
#chr1    16      +       ACCCT   0.000   1.00    0       1       8       8       0.000   0.793                                             #
#2. Must provide a chromosome_size_file                                                                                                    #
#chr1        197195432                                                                                                                          #
#chr2        181748087                                                                                                                          #
#3.
############################################################################################################################################
##########options##################
usage="usage: %prog [options], e.g. DMRcalling -m mutant_methratio.txt -c control_methratio.txt -t CG -r tair10_chrom_size.txt -b 100 -v 4 -o control_vs_mutant_DMR.txt -l 4 -d 0.4 -p 0.05 || DMRcalling is developed by Wanlu Liu from Jacobsen lab, UCLA"
parser = optparse.OptionParser(usage=usage,version="%prog 1.1")
parser.add_option("-m","--meth",dest="meth",help='methratio.txt',metavar="FILE")
parser.add_option("-r", "--bed", dest="ref", metavar="FILE", help="retion.bed (required, must be sorted)", default="")
parser.add_option("-t", "--type", dest="type", type="string", help="Cytosine context, CG/CHG/CHH", default="CG")
parser.add_option("-o", "--out", dest="outfile", metavar="FILE", help="output file title. (required)", default="mutant_vs_control_type_DMR.txt")

options, infiles = parser.parse_args()
if len(sys.argv[1:]) == 0:
    print "no argument given!"
    parser.print_help()

#if len(infiles) == 0: parser.error("Require at least methratio files for mutant (use -m) and control(use -c), and chromosome size file (-r) ")

methratio=open(options.meth,'r')
reference=open(options.ref,'r')

def cytosine_type(methratio,type):
    import re
    SearchStr= '(\S)(\S)(\S)(\S)(\S)'
    LineNumber,Result=0,0
    Num1,Num2,Num3,Num4,Num5=1,2,3,4,5
    methchr,methpos,methccount,methctcount=[],[],[],[]
    #CG
    for Line1 in methratio:
        #print(Line1)
        Line1=Line1.strip('\n')
        ElementList=Line1.split('\t')
        #print(ElementList[3])
        if (len(ElementList[3])>0):
            Result=ElementList[3]
            if type=="CG":
                if ElementList[3]== 'CG' :
                    methchr.append(ElementList[0])
                    methpos.append(ElementList[1])
                    methccount.append(ElementList[4])
                    methctcount.append(ElementList[5])
            #                   context.append(ElementList[3])
            if type=="CHG":
                if ElementList[3]== 'CHG' :
                    methchr.append(ElementList[0])
                    methpos.append(ElementList[1])
                    methccount.append(ElementList[4])
                    methctcount.append(ElementList[5])
            #                   context.append(ElementList[3])
            if type=="CHH":
                if ElementList[5]== 'CHH' :
                    methchr.append(ElementList[0])
                    methpos.append(ElementList[1])
                    methccount.append(ElementList[4])
                    methctcount.append(ElementList[5])
            #                   context.append(ElementList[3])
#context.append(ElementList[3])
        LineNumber=LineNumber+1
    return(methchr,methpos,methccount,methctcount)
def extract_methratio(methchr,methpos,methccount,methctcount,chr,start,end):
    ccount_bin,ctcount_bin=[],[]
    ## get offset of methchr
    offset = {}
    length = {}
    i = 0
    #print(i)
    current_chr = methchr[i]
    offset[current_chr] = i
    current_length = 0
    while (i < len(methchr) - 1):
        #print(i)
        if current_chr != methchr[i+1]:
            length[current_chr] = current_length
            current_length = -1
            current_chr = methchr[i+1]
            offset[current_chr] = i+1
        i = i + 1
        current_length = current_length + 1
    length[current_chr] = current_length
    for i in range(len(start)):
        current_start = start[i]
        current_end = end[i]
        current_chr = chr[i]
        if current_chr in offset:
            left = offset[current_chr] - 1
            right = offset[current_chr] + length[current_chr] + 1
            #print(left, right)
            while (left + 1 != right):
                middle = left + int((right - left) / 2)
                #print(left, middle, right)
                if (int(methpos[middle]) < int(current_start)):
                    left = middle
                #print(middle)
                else:
                    right = middle
            if right >= offset[current_chr] + length[current_chr] + 1:
                ctcount_bin.append(0)
                ccount_bin.append(0)
            #covfilter.append(0)
            else:
                j = right
                flag=0
                ccount_sum = 0
                ctcount_sum = 0
                while (int(methpos[j]) <= int(current_end)) and (j < len(methpos)) and (methchr[j] == current_chr):
                    #and (int(methctcount[j])>=int(filter)):
                    #ccount_bin.append(methcount[j])
                    #ctcount_bin.append(methccount[j])
                    ccount_sum = ccount_sum + int(methccount[j])
                    ctcount_sum = ctcount_sum + int(methctcount[j])
                    flag=flag+1
                    j+=1
                    if (j >= len(methpos)):
                        break
                ctcount_bin.append(ctcount_sum)
                ccount_bin.append(ccount_sum)
#   covfilter.append(flag)
        else:
            ctcount_bin.append(0)
            ccount_bin.append(0)
 #           covfilter.append(0)
    return(ccount_bin,ctcount_bin)#,covfilter)

outputname=options.outfile
##step1 extract methylaiton type
methchr,methpos,methccount,methctcount=cytosine_type(methratio,options.type)
print("Finish extract %s from %s"%(options.type,methratio))
step1_time=time.time()
print("Current memory usage is %s Mb" %(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024*1024)))
print("---Step 1 takes %s seconds ---" % (round(step1_time - start_time,3)))
#del(mutant,control)
#print([name for name in vars().keys()])

##step2 get range from bed file
chr,start,end=[],[],[]
for Line in reference:
    #print(Line1)
    Line=Line.strip('\n')
    ElementList=Line.split('\t')
    chr.append(ElementList[0])
    start.append(ElementList[1])
    end.append(ElementList[2])
print("---Step 2 complete")

##step4 extract methratio into bin
meth_ccount_bin,meth_ctcount_bin=extract_methratio(methchr,methpos,methccount,methctcount,chr,start,end)
print("---Step 3 complete")
#print(meth_ccount_bin)
#print(meth_ctcount_bin)
#print(len(chr),len(start),len(end),len(meth_ccount_bin),len(meth_ctcount_bin))

header="Chr\tstart\tend\tC_count\tCT_count"
numpy.savetxt(outputname, numpy.c_[chr,start,end,meth_ccount_bin,meth_ctcount_bin],fmt=('%s','%s','%s','%s','%s'),header=header,comments='',delimiter='\t')
#print("Finish %s DMR calling for %s vs %s"%(options.type,mutname,ctrlname))
#print([name for name in vars().keys()])
print("---Total time is %s seconds ---" % (round(time.time() - start_time,3)))