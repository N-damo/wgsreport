import os
import sys
import pandas as pd
from metrics_search import Search, nesteddict
from linecache import getline
from collect_insert_size import CollectInsertSizeMetrics


class MarkDuplicate(CollectInsertSizeMetrics):
    def __init__(self, module):
        self.module = module
        self.collect_file = Search('MarkDuplicate').match()
        self.create_log()

    def dedup_log(self):
        dedup = nesteddict()
        for sample in self.collect_file:
            log = self.collect_file[sample]
            # print(log)
            line = self.get_line(log, 'LIBRARY') + 1
            rec = getline(log, line).strip().split()
            UNPAIRED_READS_EXAMINED, READ_PAIRS_EXAMINED, UNMAPPED_READS, PERCENT_DUPLICATION, ESTIMATED_LIBRARY_SIZE = rec[
                1], rec[2], rec[4], rec[-2], rec[-1]
            #sample[sample]['Unmerged_reads'] = UNPAIRED_READS_EXAMINED
            #sample[sample]['Merged_reads'] = READ_PAIRS_EXAMINED
            dedup[sample]['Mapped_rate'] = round(
                (int(ESTIMATED_LIBRARY_SIZE)-int(UNMAPPED_READS))/int(ESTIMATED_LIBRARY_SIZE), 4)
            dedup[sample]['PCR_duplication'] = PERCENT_DUPLICATION
            #sample[sample]['Library_size'] = ESTIMATED_LIBRARY_SIZE
        df = pd.DataFrame.from_dict(dedup, orient='index')
        df = df.reset_index()
        df = df.rename(columns={'index': 'Samples'})
        return df


if __name__ == '__main__':
    MarkDuplicate('2.mapping').dedup_log()
