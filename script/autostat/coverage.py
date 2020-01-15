import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
from metrics_search import Search, nesteddict
from linecache import getline
from collect_insert_size import CollectInsertSizeMetrics
import numpy as np
import seaborn as sns


class CollectWgsMetrics(CollectInsertSizeMetrics):
    def __init__(self, module):
        self.module = module
        self.collect_file = Search('CollectWgsMetrics').match()
        self.create_log()

    def parse_metric(self):
        coverage = nesteddict()
        for sample in self.collect_file:
            log = self.collect_file[sample]
            line = self.get_line(log, 'GENOME_TERRITORY') + 1
            rec = getline(log, line).strip().split()
            MEAN_COVERAGE, MEDIAN_COVERAGE, PCT_1X, PCT_5X, PCT_10X, PCT_20X, PCT_30X = rec[
                1], rec[3], rec[12], rec[13], rec[14], rec[16], rec[18]
            coverage[sample]['MEAN_COVERAGE'] = MEAN_COVERAGE
            coverage[sample]['MEDIAN_COVERAGE'] = MEDIAN_COVERAGE
            coverage[sample]['Coverage_at_least_1X'] = PCT_1X  # 即测序覆盖度
            coverage[sample]['Coverage_at_least_5X'] = PCT_5X
            coverage[sample]['Coverage_at_least_10X'] = PCT_10X
            coverage[sample]['Coverage_at_least_20X'] = PCT_20X
            coverage[sample]['Coverage_at_least_30X'] = PCT_30X
            line = self.get_line(log, 'coverage')
            result = self.read_his(log, line)
            histogram = pd.DataFrame.from_dict(result, orient='index').T
            histogram['coverage'] = histogram.index
            histogram['count'] = histogram['count'].astype(int)
            histogram['coverage'] = histogram['coverage'].astype(int)
            histogram = histogram.sort_values('coverage')
            histogram.to_csv('{module}/{sample}/depth_distribution.csv'.format(
                module=self.module, sample=sample), index=False)
            self.plot_his(histogram, sample)
            self.plot_cdf(histogram, sample)
            self.genome_cov(sample)

        df = pd.DataFrame.from_dict(coverage, orient='index')
        df = df.reset_index()
        df = df.rename(columns={'index': 'Samples'})
        return df

    def plot_his(self, histogram, sample):
        histogram['count'] = histogram['count']/sum(histogram['count'])
        plt.clf()
        plt.plot('coverage', 'count', data=histogram,
                 color='mediumvioletred', linestyle='solid')
        plt.title('depth_distribution')
        plt.xlabel('depth')
        plt.ylabel('percent of reads(%)')
        plt.savefig('{module}/{sample}/depth_distribution.png'.format(
            module=self.module, sample=sample), format='png', dpi=200)
        plt.close()

    def plot_cdf(self, histogram, sample):
        histogram = histogram.sort_values('coverage')
        histogram = histogram.reset_index()
        histogram['cumulative'] = 0
        for i in range(len(histogram)):
            histogram.loc[i, 'cumulative'] = sum(histogram.loc[:i, 'count'])
        histogram['cumulative'] = histogram['cumulative'] / \
            sum(histogram['count'])
        plt.clf()
        plt.plot('coverage', 'cumulative', data=histogram,
                 color='blue', linestyle='solid')
        plt.title('depth_cumulative.png')
        plt.xlabel('coverage')
        plt.ylabel('acumulative reads proportion(%)')
        plt.savefig('{module}/{sample}/depth_cumulative.png'.format(
            module=self.module, sample=sample), format='png')
        plt.close()
        histogram = histogram[['coverage', 'cumulative']]
        histogram.to_csv('{module}/{sample}/depth_cumulative.csv'.format(
            module=self.module, sample=sample), index=False, header=True)

    def genome_cov(self, sample):
        print(sample)
        df = pd.read_csv('{}_bedgraph.bed'.format(
            sample), sep='\t', header=None)
        df.columns = ['chr', 'start', 'end', 'depth']
        df['depth'] = df['depth'].map(lambda x: np.log2(x+(1e-3)))
        plt.clf()
        plt.style.use('seaborn-darkgrid')
        num = 0
        fig = plt.figure(figsize=(10, 20), tight_layout=True)
        for chr in df['chr'].drop_duplicates():
            num += 1
            plt.subplot(12, 2, num)
            data = df[df['chr'] == chr]
            plt.plot(range(len(data)), data['depth'], label=chr,
                     marker='', color='green', linewidth=1.9, alpha=0.9)
            plt.axis(ymin=0, xmin=0)
            plt.fill_between(
                range(len(data)), data['depth'], data=df[df['chr'] == chr], color='green')
            plt.title(chr, loc='left', fontsize=12,
                      fontweight=0, color='black')
        # plt.tight_layout()
        plt.suptitle("the Distribution of Genome Coverage", fontsize=13,
                     fontweight=0, color='black', style='italic', y=1.02)
        fig.text(0.5, 0, 'Chromosome Postion(100kb)', ha='center', va='center')
        fig.text(0, 0.5, 'Average Depth_log2', ha='center',
                 va='center', rotation='vertical')
        plt.savefig('{module}/{sample}/genome_cov.png'.format(module=self.module,
                                                              sample=sample), format='png', bbox_inches='tight', dpi=200)
        plt.close()


if __name__ == '__main__':
    CollectWgsMetrics('2.mapping').parse_metric()
