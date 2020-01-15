#!/usr/bin/env python3
# coding: utf-8


import os
import sys
import pandas as pd


def get_ref(pop):
    d = []
    with open(pop, 'rt') as f:
        for i in f:
            d.append(i.strip())
    return d


def get_fam(fam):
    d = []
    with open(fam, 'rt') as f:
        for i in f:
            d.append(i.strip().split()[0])
    return d


def combine_ref_fam(pop, fam):
    pop = get_ref(pop)
    fam = get_fam(fam)
    m = zip(pop, fam)
    row = []
    for i in m:
        region, sample = i
        if region == '-':
            row.append(sample)
        else:
            row.append(region)
    return row


def cal(admix, row):
    df = pd.read_csv(admix, sep=' ', header=None)
    kilogenome=pd.read_table('/share/data3/lianlin/admixture/1kgenome.txt',sep='\t',header=None)
    panel=set(kilogenome[0])
    df.insert(0, 'Samples', row)
    df=df.drop_duplicates()
    df=df.reset_index()
    df=df.drop('index',axis=1)
    for i in range(len(df)):
        pop=df.loc[i,'Samples']
        if pop in panel:
            value=list(df.loc[i,df.columns[1:]])
            index=value.index(0.99975)
            df=df.rename(columns={index:pop})

    samples = set(df['Samples']) - set(panel)
    df2 = df[df['Samples'].isin(samples)]
    return df2


def avg_chromosome(row):
    df = cal('panel_prune.1kg.chr1.26.Q', row)
    panel = df.columns[1:]
    samples = df['Samples']
    for i in range(2, 23):
        df2 = cal('panel_prune.1kg.chr{}.26.Q'.format(i), row)
        df = df[panel]+df2[panel]
    df = df/22
    df.insert(0, 'Samples', samples)
    if os.path.exists('admixture'):
        pass
    else:
        os.mkdir('admixture')
    df.to_csv('admixture/admixture.1kg.26public.csv', index=False)


if __name__ == '__main__':
    row = combine_ref_fam('panel_prune.1kg.pop',
                          'panel_prune.1kg.fam')

    avg_chromosome(row)
    #cal('panel_prune.1kg.chr1.26.Q', row)