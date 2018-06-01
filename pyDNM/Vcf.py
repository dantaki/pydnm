#!/usr/bin/env python3
from pyDNM.Backend import err_fh,tokenize,try_index
from pyDNM.Features import Feature
import numpy as np
import sys
def genotype_pl_index(gt=None):
    """
        In presence of the GT field the same
        ploidy is expected and the canonical order is used; without GT field, diploidy is assumed.  If A is the allele in
        REF and B,C,...  are the alleles as ordered in ALT, the ordering of genotypes for the likelihoods is given by:
        
        F(j/k) = (k*(k+1)/2)+j 
    """
    _a = gt.replace('|','/').split('/')
    if len(_a)==2:
        a1,a2 = int(_a[0]),int(_a[1])
        t = a1
        if a2 < a1: a1,a2=a2,t
        return int((a2*(a2+1)/2.)+a1)
    else:
        return -1
class Vcf():
    def __init__(self):
        self.ids={} # [iid] = column index
        self.format=None # [format] = index
        self.gt=None # [iid] = genotype
        self.missing=None# True if one sample has a missing genotype
    def index_samples(self,l=None,Fam=None):
        """
        Processes CHROM header line
        stores column indices with sample ID
        """
        r = tokenize(l)
        if r==0:
            print('FATAL ERROR: unable to tokenize header {}\n'.format(l))
            sys.stderr.write('FATAL ERROR: unable to tokenize header {}\n'.format(l))
            sys.exit(1)
        for i in range(9,len(r)):
            if Fam.sex.get(r[i])==None: continue
            self.ids[r[i]]=i
    def index_format(self,l=None):
        """
        Determines indices of FORMAT values
        """
        self.format={}
        r = l.split(':')
        for i in range(0,len(r)): self.format[str(r[i])]=i
    def load_genotypes(self,r=None):
        """
        Stores genotypes
        """
        self.missing = [False, []]
        self.gt={}
        for iid in self.ids:
            self.gt[iid]=str(r[self.ids[iid]].split(':').pop(self.format['GT']))
        for x in self.gt:
            if '.' in self.gt[x]:
                self.missing[0]=True
                self.missing[1].append('{}:{}'.format(x,self.gt[x]))
    def check_genotypes(self,iid=None,variant=None):
        """
        Logic to check for genotype value in dict
        """
        gt=None
        if self.gt.get(iid)!=None: gt=self.gt[iid]
        if gt==None:
            sys.stderr.write('WARNING: missing genotype entry: {} {}\n'.format(variant,iid))
        return gt
    def allele_depth(self,entry=None,dnm=None,par=None):
        """
        returns 
            the allele ratio of the dnm and parent alleles
            the log2 coverage ratio (relative to median)
        """
        _buff = 1 # add one to avoid 0 values: taken from forestDNM
        _ad = str(entry).split(':')
        if len(_ad)<self.format['AD']: return -1,np.nan
        else:
            _ad = _ad[self.format['AD']].split(',')
            _dnm = try_index(_ad,dnm)
            _par = try_index(_ad,par)
            if _dnm==None or _par==None: return -1,np.nan
            else: 
                _par = float(_par)+_buff
                ar = float(_dnm)/ (_par)
                med = np.median([float(x) for x in _ad if x!='.']) 
                if not np.isfinite(med): return -1, np.nan
                else: return ar, np.log2( (float(_dnm)+float(_par)) / (med+_buff) )
    def genotype_quals(self,entry=None):
        """
        returns GQ values
        """
        _gq = str(entry).split(':')
        if len(_gq)<self.format['GQ']: return -1
        else:
            _gq = _gq[self.format['GQ']]
            if _gq == '.': return -1
            else: return float(_gq) 
    def phred_quals(self,entry=None,gt=None):
        """
        returns PL values
        """
        _pl = str(entry).split(':')
        if len(_pl)<self.format['PL']: return -1
        else:
            _pl = _pl[self.format['PL']].split(',')
            idx = genotype_pl_index(gt)
            if idx==-1: return -1
            else:
                pl = try_index(_pl,idx)
                if pl==None: return -1
                elif pl=='.': return -1
                else: return float(pl)
    def parse(self,fh=None, Fam=None,verb=None,ofh=None):
        err_fh(fh)
        vcf_fh=None
        if fh.endswith('.gz'):
            import gzip
            vcf_fh = gzip.open(fh,'rt', 9)
        else: vcf_fh = open(fh,'r')
        if vcf_fh==None:
            print('FATAL ERROR: {} file format unknown.\nAccepted VCF formats: text, bgzip (ends with \'.gz\')\n'.format(fh))
            sys.stderr.write('FATAL ERROR: {} file format unknown.\nAccepted VCF formats: text, bgzip (ends with \'.gz\')\n'.format(fh))
            sys.exit(1)
        out = open(ofh,'w')
        out.write('chrom\tpos\tid\tref\talt\tiid\toffspring_gt\tfather_gt\tmother_gt\t')
        out.write('{}\n'.format(Feature().header()))
        for l in vcf_fh:
            if l.startswith('#CHROM'): self.index_samples(l,Fam)
            if l.startswith('#'): continue
            r=tokenize(l)
            if r==0: continue
            variant = tuple(map(str,r[0:5]))
            # load genotypes
            self.index_format(r[8])
            self.load_genotypes(r)
            if self.missing[0]==True:
                if verb==True: sys.stderr.write('WARNING: missing genotypes {} {}\n'.format(variant,','.join(self.missing[1])))
                continue
            """
            Foreach trio
            """
            for kid in Fam.offspring:
                if self.ids.get(kid)==None: continue
                dad,mom = Fam.offspring[kid]
                kgt,dgt,mgt = self.check_genotypes(kid,variant),self.check_genotypes(dad,variant),self.check_genotypes(mom,variant)
                # skip if genotypes are not available
                if kgt==None or dgt==None or mgt==None: continue
                alleles = kgt.replace('|','/').split('/')
                """
                determine which of the offspring alleles are de novo and inherited
                """
                par,dnm = None,None # parent allele, de novo allele
                for a in alleles:
                    # if the offspring allele is NOT in parents genotype ==> De Novo
                    if str(a) not in dgt and str(a) not in mgt: dnm=int(a)
                    # if the offspring allele is in ONE of the parents genotype ==> Inherited
                    elif str(a) in dgt or str(a) in mgt: par=int(a)

                # skip if there is not one inherited and one de novo allele
                if par==None or dnm==None: continue
                
                # init Features
                Feat = Feature()
                # load INFO features
                Feat.parse(r)
                # Allele depth
                if self.format.get('AD')==None: 
                    if verb==True: sys.stderr.write('WARNING: missing allele depth AD {}\n'.format(variant))
                    continue
                dad_ad,dad_dp = self.allele_depth(r[self.ids[dad]],dnm,par)
                mom_ad,mom_dp = self.allele_depth(r[self.ids[mom]],dnm,par)
                kid_ad,kid_dp = self.allele_depth(r[self.ids[kid]],dnm,par)
                if dad_ad==-1 or mom_ad==-1 or kid_ad==-1: continue
                if not np.isfinite(dad_dp) or not np.isfinite(mom_dp) or not np.isfinite(kid_dp): continue
                Feat.p_ar_max= max(dad_ad,mom_ad)
                Feat.p_ar_min= min(dad_ad,mom_ad)
                Feat.o_ar = kid_ad
                Feat.p_dp_max = max(dad_dp,mom_dp)
                Feat.p_dp_min = min(dad_dp,mom_dp)
                Feat.o_dp = kid_dp
                # Genotype quality scores
                if self.format.get('GQ')==None:
                    if verb==True: sys.stderr.write('WARNING: missing genotype quality GQ {}\n'.format(variant))
                    continue
                kid_gq,dad_gq,mom_gq = self.genotype_quals(r[self.ids[kid]]),self.genotype_quals(r[self.ids[dad]]),self.genotype_quals(r[self.ids[mom]])
                if kid_gq==-1 or dad_gq==-1 or mom_gq==-1: continue
                Feat.p_gq_max = max(dad_gq,mom_gq)
                Feat.p_gq_min = min(dad_gq,mom_gq)
                Feat.o_gq = kid_gq
                Feat.n_alt = len(str(r[4]).split(','))
                # Phred scaled genotype likelihoods
                if self.format.get('PL')==None:
                    if verb==True: sys.stderr.write('WARNING: missing Phred-adjusted genotype likelihoods PL {}\n'.format(variant))
                    continue
                kid_pl, dad_pl, mom_pl = self.phred_quals(r[self.ids[kid]],kgt),self.phred_quals(r[self.ids[dad]],dgt),self.phred_quals(r[self.ids[mom]],mgt)
                if kid_pl==-1 or dad_pl==-1 or mom_pl==-1: continue
                kid_d_pl,kid_m_pl, dad_o_pl, mom_o_pl = self.phred_quals(r[self.ids[kid]],dgt),self.phred_quals(r[self.ids[kid]],mgt),self.phred_quals(r[self.ids[dad]],kgt),self.phred_quals(r[self.ids[mom]],kgt)
                if kid_d_pl==-1 or kid_m_pl==-1 or dad_o_pl==-1 or mom_o_pl == -1: continue
                Feat.p_og_max = max(dad_o_pl,mom_o_pl)
                Feat.p_og_min = min(dad_o_pl,mom_o_pl)
                Feat.p_pg_max = max(dad_pl,mom_pl)
                Feat.p_pg_min = min(dad_pl,mom_pl)
                Feat.og = kid_pl
                Feat.o_pg = np.median([kid_d_pl,kid_m_pl])
                o = Feat.output()
                out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format('\t'.join(variant),kid,kgt,dgt,mgt,o))
        out.close()