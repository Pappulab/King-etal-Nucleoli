#!/usr/bin/env python
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pprint import pprint
import scipy.stats as st
from tabulate import tabulate


sns.set_style('whitegrid')


# Supremely useful answer to convert a DataFrame to a TSV file: https://stackoverflow.com/a/35974742/866930
def to_fwf(df, fname):
    content = tabulate(df.values.tolist(), list(df.columns), tablefmt="plain", floatfmt='.7f')
    open(fname, 'w').write(content)
pd.DataFrame.to_fwf = to_fwf


def main():
    confidence_interval = 0.95
    directory = os.path.dirname(os.path.abspath(__file__))
    files = [os.path.join(directory, p) for p in sorted(os.listdir(directory)) if p.endswith('.txt')]
    condensates = list(sorted(f for f in files if 'condensate' in f and 'crop' not in f.lower()))
    patches = list(sorted(f for f in files if 'patch' in f and 'crop' not in f.lower()))
    for c, p in zip(condensates, patches):
        cs = os.path.basename(c.replace('_condensates.txt', ''))
        ps = os.path.basename(p.replace('_patches.txt', ''))
        assert cs == ps

    classes = list(set(['_'.join(os.path.basename(f).replace('_condensates.txt', '').split('_')[:-1]) for f in condensates if 'crop' not in f.lower()]))

    condensates_medians = list()
    patches_medians = list()
    data = dict()
    for class_name in classes:
        c_medians_name = list()
        p_medians_name = list()
        for c, p in zip(condensates, patches):
            c_medians = np.loadtxt(c, skiprows=1)

            pf = pd.read_csv(p, sep='\s+')
            p_medians = list()
            for name in pf.columns.values:
                p_median = np.median(pf[name])
                p_medians.append(p_median)

            if class_name in c and class_name in p:
                c_medians_name += c_medians.tolist()
                p_medians_name += p_medians
        _dilute = np.array(p_medians_name).flatten()
        _condensate = np.array(c_medians_name).flatten()
        data[class_name] = dict()
        data[class_name]['condensate'] = _condensate[~np.isnan(_condensate)]
        data[class_name]['dilute'] = _dilute[~np.isnan(_dilute)]

        for key in data[class_name]:
            values = data[class_name][key].tolist()
            # Choose the sampling method to calculate the confidence interval based on the number of points
            if len(values) < 30:
                confidence = st.t.interval(alpha=confidence_interval, df=len(values)-1, loc=np.mean(values), scale=st.sem(values))
            else:
                confidence = st.norm.interval(alpha=confidence_interval, loc=np.mean(values), scale=st.sem(values))

            # points within confidence interval
            confident_values = [v for v in values if confidence[0] <= v <= confidence[1]]
            data[class_name][key] = confident_values

    colors = 'blue red'.split()

    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_subplot(111)
    all_dfs = list()
    hue = list()
    for class_name in data:
        d = dict()
        d[f'{class_name}-condensate'] = pd.Series(data[class_name]['condensate'])
        d[f'{class_name}-dilute'] = pd.Series(data[class_name]['dilute'])
        hue.append('condensate')
        hue.append('dilute')

        df = pd.DataFrame(d)
        all_dfs.append(df)
    df = pd.concat(all_dfs)
    df.to_fwf('all_data.tsv')


    ndf = dict()
    ndf['Name'] = list()
    ndf['Dense'] = list()
    ndf['Dilute'] = list()
    for class_name in data:
        fig = plt.figure(figsize=(6, 8))
        ax = fig.add_subplot(111)

        d = dict()
        d['Dense'] = pd.Series(data[class_name]['condensate'])
        d['Dilute'] = pd.Series(data[class_name]['dilute'])

        ndf['Name'].append(class_name)
        ndf['Dense'].append(d['Dense'].median())
        ndf['Dilute'].append(d['Dilute'].median())

        colors = 'blue red blue red blue red blue red'.split()
        g = sns.stripplot(df, ax=ax, palette=colors, zorder=-1)
        g = sns.violinplot(df, ax=ax, palette=colors, zorder=1)
        for points in g.collections[::1]:
            points.set_alpha(0.6)
        g.set_xticklabels(g.get_xticklabels(), rotation=45)

        ax.set_xlabel('Measurement')
        ax.set_ylabel('Calculated pH')
        [x.set_linewidth(3) for x in ax.spines.values()]
        fig.tight_layout()
        plt.show()
    pprint(ndf)



if __name__ == '__main__':
    main()

