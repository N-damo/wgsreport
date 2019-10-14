#!/usr/bin/env python2
#coding:utf-8

import json
from collections import defaultdict

def f():
    pop_dic=defaultdict(list)
    with open('/share/data3/lianlin/soft/bin/wgs/1kgenome.txt','r') as f:
        for i in f:
            line=i.strip().split()
            pop=line[0]
            sample=line[1]
            site=line[2]
            pop_dic[site].append(sample)
    with open('1kgenome.json','w') as f:
        json.dump(pop_dic,f)
f()
