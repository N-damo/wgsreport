import os
import sys



def merge_fam(fam):
    d=[]
    with open(fam,'rt') as f:
        for k,v in enumerate(f):
            line=v.strip().split()
            sample=line[0]
            d.append(sample)
    return d

def sample_map(ref):
    d={}
    with open(ref,'rt') as f:
        for k,v in enumerate(f):
            line=v.strip().split()
            sample=line[1]
            country=line[0]
            d[sample]=country
    return d


def pop(fam,ref,prefix):
    fam=merge_fam(fam)
    ref=sample_map(ref)
    with open(prefix+'.pop','wt') as f:
        for sample in fam:
            country=ref.get(sample,'-')
            #f.write(sample+'\t'+country+'\n')
            f.write(country+'\n')


if __name__ == '__main__':
    ref=sys.argv[1]
    fam=sys.argv[2]
    prefix=sys.argv[3]
    pop(fam,ref,prefix)





