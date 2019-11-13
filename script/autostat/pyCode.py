#!/usr/bin/env python3
#coding:utf-8

import sys 
import os
import subprocess
from linecache import getline
"""
Nozzle.R1
<meta http-equiv=”Content-Type” content=”text/html; charset=utf-8″>
"""

def convert_file_to_utf8(filename):
    utf="""<meta http-equiv=”Content-Type” content=”text/html; charset=utf-8″>"""
    with open(filename,'rt') as fin:
        with open('utf8_'+filename,'wt') as fout:
            for key,value in enumerate(fin):
                if key == 2:
                    if value.find('utf-8') != -1:
                        raise SystemError('exists\n')
                    else:
                        print('add coding utf-8')
                        fout.write(utf+'\n')
                        fout.write(value)
                else:
                    fout.write(value)


if __name__ == '__main__':
    convert_file_to_utf8(sys.argv[1])
