#!/usr/bin/env python3
# coding:utf-8

import os
import sys
import subprocess
from collections import defaultdict
from pipeline.shortVariant import ShortVariant
from pipeline.cnv import CNV
from pipeline.sv import SV
import argparse
from argparse import ArgumentParser
import logging


def get_args():
    parser = ArgumentParser(
        description='wgs workflow,including pipeline,autostat and html report generation')
    parser.add_argument('--input', '-i', help='sample info',
                        type=extend_file, required=True)
    parser.add_argument('--output_path', '-o', help='output directory',
                        default='./')
    parser.add_argument('--run', '-r', choices={
                        'on', 'off'}, help='if you choose run on,please add -r argument.the default is off', default='off')
    args = parser.parse_args()
    return args


def extend_file(file):
    if os.path.exists(file):
        pass
    else:
        raise argparse.ArgumentTypeError(
            '{file} dose not exist'.format(file=file))
    return file


class WGS(object):

    def __init__(self, input):
        self.input = input
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
                        logging.debug('it seem has less or more than five items in the {} row,it may be a blank row,so pass it and continue\n'.format(
                            row_line))
                else:
                    self.working_space = i.strip().split()
                    assert len(
                        self.working_space) == 1, 'please use batch number,e.g. gmb1,gmb2 ...'
                    self.working_space = self.working_space[0]
                    self.create_working_space()
                    self.working_space = os.path.abspath(self.working_space)
                    if os.path.exists(os.path.join(self.working_space, 'workspace')):
                        pass
                    else:
                        os.mkdir(os.path.join(self.working_space, 'workspace'))
                    row_line += 1

    def create_working_space(self):
        if not os.path.exists(self.working_space):
            logging.debug('create {} directory'.format(self.working_space))
            os.mkdir(self.working_space)
        else:
            pass


if __name__ == '__main__':
    args = get_args()
    input = args.input
    run = args.run
    output_path=args.output_path#project 目录
    logging.basicConfig(level=logging.DEBUG,filename=os.path.join(output_path,'app.log'),
                        format='%(levelname)s:%(asctime)s:%(message)s')
    os.chdir(os.path.dirname(output_path))
    pipeline = WGS(input)
    # 样品信息，key=sampleName,value='fq1,fq2,adapter,library,sampleName
    groups = pipeline.groups
    working_space = pipeline.working_space  # first line in input
    if run == 'off':
        logging.debug(
            'you choose dry run,so only short variant bat file will create in {}'.format(working_space))
        ShortVariant(groups, working_space).main()
        CNV(groups, working_space).CnvnatorCNV()
        SV(groups, working_space).dellySV()
    else:
        ShortVariant(groups, working_space).main('on')
        CNV(groups, working_space).CnvnatorCNV(run='on')
        SV(groups, working_space).dellySV(run='on')
