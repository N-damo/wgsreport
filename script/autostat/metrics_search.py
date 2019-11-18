#!/usr/bin/env python3
# coding:utf-8

import os
import sys
from collections import defaultdict


class Search(object):

    module_key = {'fastp': 'json', 'collectinsertsizemetrics': 'insert_size_metrics.txt',
                  'markduplicate': 'dedup.log', 'collectwgsmetrics': 'wgs_metrics.txt'}

    def __init__(self, module):
        self.module = module.lower()
        if self.module in self.module_key:
            pass
        else:
            raise SystemError('{} not in matched module list.\n I support module list is {}.\n '.format(
                self.module, list(self.module_key.keys())))

    def match(self):
        collect_file = {}
        for directory, _, file_list in os.walk('./'):
            for file in file_list:
                if file.endswith(self.module_key[self.module]):
                    sample_name = file.split('.')[0]
                    collect_file[sample_name] = os.path.abspath(
                        os.path.join(directory, file))
        if len(collect_file) == 0:
            raise SystemError(
                'can not find {} log file.\n'.format(self.module))
        else:
            pass
        return collect_file


def nesteddict():
    return defaultdict(nesteddict)
