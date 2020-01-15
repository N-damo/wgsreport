import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
from metrics_search import Search, nesteddict
from linecache import getline


class CollectInsertSizeMetrics(object):
    def __init__(self, module):
        self.module = module
        # +Search('MarkDuplicate').match()+Search('CollectWgsMetrics').match()
        self.collect_file = Search('CollectInsertSizeMetrics').match()
        self.create_log()
        # self.insert_size()

    def create_log(self):
        if os.path.exists(self.module):
            pass
        else:
            os.mkdir(self.module)
        for sample in self.collect_file:
            direc = os.path.join(self.module, sample)
            if os.path.exists(direc):
                pass
            else:
                os.mkdir(direc)

    def get_line(self, file, pattern):
        with open(file, 'rt') as f:
            for index, value in enumerate(f):
                if value.startswith(pattern):
                    return index + 1  # 返回下标1开始的
                else:
                    continue

    def read_his(self, file, index):
        result = nesteddict()
        with open(file, 'rt') as f:
            for key, value in enumerate(f):
                if key < index:
                    continue
                else:
                    try:
                        cov, count = value.strip().split()
                        result['count'][cov] = count
                    except ValueError:
                        break

        return result

    def insert_size(self):
        insert = nesteddict()
        for sample in self.collect_file:
            log = self.collect_file[sample]
            line = self.get_line(log, 'MEDIAN_INSERT_SIZE') + 1  # 从1开始
            rec = getline(log, line).strip().split()
            insert_mean, insert_std = rec[5], rec[6]
            insert[sample]['InsertSize_mean'] = insert_mean
            insert[sample]['Insert_std'] = insert_std
            line = self.get_line(log, 'insert_size')
            result = self.read_his(log, line)
            histogram = pd.DataFrame.from_dict(result, orient='index').T
            histogram['insert_size'] = histogram.index
            # histogram=histogram.rename(columns={'count':'Forward_count'})
            histogram['count'] = histogram['count'].astype(int)
            histogram['insert_size'] = histogram['insert_size'].astype(int)
            histogram = histogram.sort_values('insert_size')
            histogram['count'] = histogram['count']/sum(histogram['count'])
            histogram.to_csv('{module}/{sample}/insert_size.csv'.format(
                module=self.module, sample=sample), index=False, header=True)
            self.plot_insert_size(histogram, sample)

        df = pd.DataFrame.from_dict(insert, orient='index')
        df = df.reset_index()
        df = df.rename(columns={'index': 'Samples'})
        return df

    def plot_insert_size(self, file, sample):
        plt.clf()
        plt.plot('insert_size', 'count', data=file,
                 color='black', linestyle='solid')
        plt.title('insert size distribution')
        plt.xlabel('insert size')
        plt.ylabel('read percent(%)')
        plt.tight_layout()
        plt.savefig('{module}/{sample}/insert_size.png'.format(
            module=self.module, sample=sample), format='png', dpi=200)
        plt.close()


if __name__ == '__main__':
    CollectInsertSizeMetrics('2.mapping').insert_size()
