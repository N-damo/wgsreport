#!/usr/bin/env python3
#coding:utf-8

import os
import sys
import subprocess
from collections import defaultdict
#from pipeline import cnv,sv,shortVariant
from pipeline.shortVariant import ShortVariant
from pipeline.cnv import CNV
from pipeline.sv import SV
from autostat import wgshsmetrics
import argparse
from argparse import ArgumentParser
from report.report import html_report


def get_args():
    parser=ArgumentParser(description='wgs workflow,including pipeline,autostat and html report generation')
    parser.add_argument('--input','-i',help='sample info',type=extend_file,required=True)
    parser.add_argment('--run','-r',choice={'on','off'},help='if you choose run on,please add -r argument.the default is off',default='off')


def extend_file(file):
    if os.path.exists(file):
        pass
    else:
        raise argparse.ArgumentTypeError('{file} dose not exist'.format(file=file))
    return file



class WGS(object):
    
    def __init__(self,input):
        self.input=input
        self.group()


    def group(self):
        self.groups = defaultdict(list)
        self.working_space = ''
        row_line = 0
        with open(self.input, 'r') as f:
            for i in f:
                if row_line > 0:
                    row_line += 1
                    line = i.strip().split()
                    if len(line) == 5:
                        sample = line[4]
                        self.groups[sample].append(line)
                    else:
                        sys.stdout.write('it seem has less or more than five items in the {} row,it may be a blank row,so pass it and continue\n'.format(
                            row_line))
                else:
                    self.working_space = i.strip().split()
                    assert len(
                        self.working_space) == 1, 'please use batch number,e.g. gmb1,gmb2 ...'
                    self.working_space = self.working_space[0]
                    self.create_working_space()
                    self.working_space = os.path.abspath(self.working_space)
                    row_line += 1



    def create_working_space(self):
        if not os.path.exists(self.working_space):
            os.mkdir(self.working_space)
        else:
            pass




if __name__ == '__main__':
    args=get_args()
    input=args.input
    run=args.run
    pipeline=WGS(input)
    groups=pipeline.groups#样品信息，key=sampleName,value='fq1,fq2,adapter,library,sampleName
    working_space=pipeline.working_space#first line in input
    if run == 'off':
        ShortVariant(groups,working_space).main('off')
    else:
        ShortVariant(groups,working_space).main('on')
        CNV(groups,working_space).CnvnatorCNV()
        SV(groups,working_space).dellySV()
        wgshsmetrics.stat(list(groups.keys()))
        html_report()
    