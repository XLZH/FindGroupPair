#!/usr/bin/env python

import os
import re
import time
import sys
import pysam
import os.path
from optparse import OptionParser

# global regular to match (eg: "\ATGTTG_" and "_CCGG") in seqname
Barcodeseq = re.compile('^\|[ATCG]+_')
Barcodefix = re.compile('_[ATCG]+/')


class Group(object):
    def __init__(self, gstart, gend, strand='*', readcount=0):
        '''reads: The sequence located in this group!
           readcount: The count of reads located in this group!
           gstart/gend: The group start and end!
        '''
        self.gstart = gstart # The Leftmost position of the Group
        self.gend = gend  # The Rightmost position of the Group
        self.readcount = readcount  # Calculate the count of all reads ('CCGG' and 'TAAT')
        self.strand = strand  # Use the strand of longseq which has fixseq of 'CCGG'
        self.reads = []  # Now only keep the reads with fixseq of 'CCGG'

    def addreads(self, read):
        '''This is used to add the reads to this group list!
        '''
        self.reads.append(read)

class Pairgroup(object):
    def __init__(self, Chr, Tag1, group1, Tag2, group2):
        ''' To sotre the information of paired group! 
        '''
        self.Chr = Chr
        self.distance = str(abs(group2.gstart - group1.gend))
        self.P1 = (Tag1, Chr, str(group1.gstart), str(group1.gend), str(group1.readcount))
        self.P2 = (Tag2, Chr, str(group2.gstart), str(group2.gend), str(group2.readcount))

class Singlegroup(object):
    def __init__(self, Chr, Tag, group):
        ''' To sotre the information of Single group! 
        '''
        self.Chr = Chr
        self.P1 = (Tag, Chr, str(group.gstart), str(group.gend), str(group.readcount))


def Writeout(result, templatefilepath, outfilepath):
    ''' Write out the bamfile(distributed by Barcode) with the filename of Bacode!
    '''
    resultlen = len(result)

    if(resultlen == 0):
        return 0
    try:
        templatefile = pysam.AlignmentFile(templatefilepath, 'rb')
        outfile = pysam.AlignmentFile(outfilepath, 'wb', template = templatefile)
    except IOError:
        print 'error: Failed to open templatefile or outfile![Writeout]'

    print 'Writeout the %s [%d]......\n' %(os.path.basename(outfilepath), resultlen)
    for element in result:
        outfile.write(element)

    templatefile.close()
    outfile.close()


def ReadsDistribute(bamfilepath, barcodepath, outfilepath):
    ''' eg: Barcode = {'BC002':TCATGTCC, 'BC003':TGTGAGGT, ...}
            Result = {'BC002':[reads1,reads2,...], 'BC003':[reads1,reads2,...], ...}
    '''
    Barcode = {}
    Result = {}
    discardreadscount = 0
    keepreadscount = 0

    try:
        bamfile = pysam.AlignmentFile(bamfilepath, 'rb')
        barcode = open(barcodepath, 'rb')
    except IOError:
        print 'error: Failed to open the input bamfile or barcodefile![ReadsDistribute]'
        sys.exit(-1)

    for line in barcode.readlines():
        tag = line.split('\t')[0]
        seq = line.split('\t')[1].rstrip()
        Barcode[seq] = tag
        Result[tag] = []

    print '%s' %(time.ctime())
    print 'Start to build the extracted list....[%s]\n' %(time.ctime())
    for read in bamfile.fetch(until_eof = True):
        if(read.is_proper_pair):
            try:
                barseq = Barcodeseq.search(read.query_name).group()[1:-1]
                bartag = Barcode[barseq]
                Result[bartag].append(read)
                keepreadscount += 1
            except:
                discardreadscount += 1
                continue
        else:
            discardreadscount += 1
            continue
    print 'The outresult list has been builded! [%s]' %(time.ctime())
    print 'Discardreadscount: %d\tKeepreadscount: %d\n' %(discardreadscount, keepreadscount)

    bamfile.close()
    barcode.close()

    if(not os.path.exists(outfilepath)):
        os.makedirs(outfilepath)

    for key in Result.keys():
        filename = key + '.bam'
        writeoutfile = os.path.join(outfilepath,filename)
        Writeout(Result[key], bamfilepath, writeoutfile)


def Extractbychr(inbamfile):
    ''' Chromosome = [('chr1',[reads1,reads2,...]),('chr2',[reads1,reads2,...]),...,('chrY',[reads1,reads2,...])]
    '''
    Chromosome = [('chr1',[])]
    LastReadChrindex = 0

    try:
        bamFile = pysam.AlignmentFile(inbamfile, 'rb')
    except IOError:
        print 'error: Failed to open %s![Extractbychr]' %(os.path.basename(inbamfile))
        sys.exit(1)

    for read in bamFile.fetch(until_eof = True):
        if(read.rname >=0 and read.rname <24):
            Chr = bamFile.getrname(read.rname)
        else:
            continue

        if(Chr == Chromosome[LastReadChrindex][0]):
            Chromosome[LastReadChrindex][1].append(read)
        else:
            Chromosome.append((Chr, [read]))
            LastReadChrindex += 1

    bamFile.close()
    return Chromosome


def Calgroupinfo(Gobject, read):
    ''' Gobject is a struct with members of gstart,gend,readcount,strand and a list of reads(with fixseq 'CCGG')!
        positiveflag/negativeflag : The flag whoes reads maped with a longer sequence (with a fixseq 'CCGG')!
    '''
    start = read.reference_start
    end = read.reference_end
    positiveflag = [99, 97]
    negativeflag = [83, 81, 89, 153]

    Gobject.gend = end if(end > Gobject.gend) else Gobject.gend
    Gobject.readcount += 1
   
    try: 
        Fixseq = Barcodefix.search(read.query_name).group()[1:-1]
        if(Fixseq == 'CCGG'):
            Gobject.addreads(read)        

            if(read.flag in positiveflag):
                Gobject.strand = '+'
            elif(read.flag in negativeflag):
                Gobject.strand = '-'
            else:
                pass
    except AttributeError:
        pass


def Classgroup(chrtuple, readdistance):
    ''' chrtuple = ('chr1', [reads1, reads2, ...])
        Groups = [('chr1','G1',group1), ..., ('chr1','G499',group499)]
        NOTE: The Groups returned only include one chromosome group information!
    '''
    Groups = []
    poshash = []
    count = -1

    for element in chrtuple[1]:
        readpos = element.pos
        grouped = 0
        i = readdistance

        if(readpos in poshash):
            poshash = range(element.reference_start, element.reference_end)
            Calgroupinfo(Groups[count][2], element)
            continue

        while(i > 0):
            testpos = readpos - i
            if(testpos in poshash):
                poshash = range(element.reference_start, element.reference_end)
                Calgroupinfo(Groups[count][2], element)
                grouped = 1
                break
            i -= 1

        if(not grouped):
            poshash = []
            count += 1
            poshash = range(element.reference_start, element.reference_end)
            groupMarker = 'G' + str(count + 1)
            Groups.append((chrtuple[0], groupMarker, Group(element.reference_start, element.reference_end)))
            Calgroupinfo(Groups[count][2], element)

    return Groups


def WriteGroupinfo(ChrsList, outgroupname, readdistance):
    ''' ChrsList = [('chr1',[reads1,reads2,...]), ('chr2',[reads1,reads2...]) ...]
    '''
    Groupinfo = []

    for chrtuple in ChrsList:
        ChrGroup = Classgroup(chrtuple, readdistance)
        for gtuple in ChrGroup:

            # struct: ginfo = ['chromosome', 'GroupTag', 'gstart', 'gend', 'readcount', 'strand']
            ginfo = [gtuple[1], gtuple[0], str(gtuple[2].gstart), str(gtuple[2].gend), str(gtuple[2].readcount), gtuple[2].strand]
            Groupinfo.append('\t'.join(ginfo))

    with open(outgroupname, 'a+') as outwritefile:
        outwritefile.write('Chr\tGroupTag\tGroupStart\tGroupEnd\tGroupReadCount\tGroupStrand\n')
        for e in Groupinfo:
            outwritefile.write('%s\n' %e)

    outwritefile.close()


def GroupFilter(ChrsList, readdistance):
    ''' Filtergroups struct as follow:
        Filtergroups = [[('chr1','G1',group1), ('chr1','G2',group2),...], ...., [('chrY','G1',group1), ('chrY','G2',group2),....]]
        NOTE:The Filtergroups contain all the group information from One barcode bam file
    '''
    Filtergroups = []

    for chrtuple in ChrsList:
        propergroup = []

        ChrGroup = Classgroup(chrtuple, readdistance)
        for gtuple in ChrGroup:
            if(gtuple[2].readcount < 4 or gtuple[2].strand == '*'):
                continue
            else:
                propergroup.append(gtuple)
        
        Filtergroups.append(propergroup)

    return Filtergroups


def CallPairinfo(chrgrouplist):
    ''' chrgrouplist = [('chr1','G1',<group1>), ('chr1','G2',<group2>), ...]
        pairsingleinfo[0] : paired groups information
        pairsingleinfo[1] : single groups information
        eg: pairsingleinfo[0] = [Pair1, Pair2, ...]
            pairsingleinfo[1] = [Single1, Single2, ...]
    '''
    pairsingleinfo = [[], []]
    marker = []

    for i in xrange(len(chrgrouplist)):
        grouped = False

        if(chrgrouplist[i][1] in marker):
            continue

        currentgroup = chrgrouplist[i][2]
        for j in xrange(i+1, len(chrgrouplist)):
            if(chrgrouplist[j][1] in marker):
                continue

            nextgroup = chrgrouplist[j][2]
            if(currentgroup.strand != nextgroup.strand):
                distance = nextgroup.gstart - currentgroup.gend

                if(distance > 4000 and distance < 15000):

                    # eg: Pairgroup('chr1', 'G1', <group1>, 'G2', <group2>)
                    pairsingleinfo[0].append(Pairgroup(chrgrouplist[i][0], chrgrouplist[i][1], chrgrouplist[i][2], chrgrouplist[j][1], chrgrouplist[j][2]))
                    marker.append(chrgrouplist[i][1])
                    marker.append(chrgrouplist[j][1])

                    grouped = True
                    break

        if(not grouped):
            pairsingleinfo[1].append(Singlegroup(chrgrouplist[i][0], chrgrouplist[i][1], chrgrouplist[i][2]))
            
    return pairsingleinfo


def WritePairinfo(pairsingleinfo, outpairname, outsinglename):
    with open(outpairname, 'a+') as pairfile:
        for pgroup in pairsingleinfo[0]:
            pinfo = '\t'.join(pgroup.P1) +'\t'+ '\t'.join(pgroup.P2) +'\t'+ pgroup.distance
            pairfile.write('%s\n' %pinfo)

    with open(outsinglename, 'a+') as singlefile:
        for sgroup in pairsingleinfo[1]:
            sinfo = '\t'.join(sgroup.P1)
            singlefile.write('%s\n' %sinfo)

    pairfile.close()
    singlefile.close()


def FindGroupPair(Filtergroups, outpairname, outsinglename):
    ''' pairsingleinfo = [[Pair1, Pair2, Pair3,...], [Single1, Single2, Single3, ...]]
    '''
    with open(outpairname, 'a+') as pairfile:
        pairfile.write('Tag1\tChr\tGstart1\tGend1\tReadCount1\tTag2\tChr\tGstart2\tGend2\tReadCount2\tDistance\n')

    with open(outsinglename, 'a+') as singlefile:
        singlefile.write('Tag\tChr\tGstart\tGend\tReadCount\n')

    pairfile.close()
    singlefile.close()

    for chrgrouplist in Filtergroups:
        pairsingleinfo = CallPairinfo(chrgrouplist)
        WritePairinfo(pairsingleinfo, outpairname, outsinglename)


def main():
    usage = 'usage: %prog [options] <input.bam> -b <barcodefile> -d <readdistacne> -o <outputdir>'
    parser = OptionParser(usage=usage)
    parser.add_option('-b', '--barcode-file', dest='barcodefile',
                        help='The barcodefile!')
    parser.add_option('-d', '--read-distance', dest='readdistance',
                        help='The distance between two reads!')
    parser.add_option('-o', '--output-dir', dest='outputpath',
                        help='The outputdir where to contain the extracted bamfile!')

    (options, args) = parser.parse_args()
    if(len(args) != 1):
        parser.print_help()
        sys.exit(0)

    inputbamFile = args[0]


    # Distribute the reads by Barcode and quality(is_proper: with tag 163, 99, 147, 83)
    #ReadsDistribute(inputbamFile, options.barcodefile, options.outputpath)


    # Class the group and write out the count file
    properbamfile = os.listdir(options.outputpath)
    Filtergruops = []

    for m in properbamfile:
        name = m.split('.')[0] + '.group'
        outgroupname = os.path.join(options.outputpath, name)
        inbamfile = os.path.join(options.outputpath, m)

        print 'Start to deal with %s! [%s]' %(os.path.basename(inbamfile), time.ctime())
        Chrs = Extractbychr(inbamfile)

        print 'Start to write the file to %s! [%s]' %(os.path.basename(outgroupname), time.ctime())
        WriteGroupinfo(Chrs, outgroupname, int(options.readdistance))

        print 'Start to Filter the Groups by readcount! [%s]' %(time.ctime())
        Filtergroups = GroupFilter(Chrs, int(options.readdistance))

        print 'Start to FinGroupinfo [%s]\n' %(time.ctime())
        pairname = m.split('.')[0] + '.pair'
        singlename = m.split('.')[0] + '.single'
        outpairname = os.path.join(options.outputpath, pairname)
        outsinglename = os.path.join(options.outputpath, singlename)

        FindGroupPair(Filtergroups, outpairname, outsinglename)


if __name__ == '__main__':
    main()
