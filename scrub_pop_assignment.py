#!/usr/bin/env python
from glob import iglob
from os import getcwd
import numpy as np
import matplotlib.pyplot as plt
import brewer2mpl

def struct_to_fig(data):
    for a in data:
        K = int(data[a]['assn'].shape[0])
        bmap = brewer2mpl.get_map('Set3', 'qualitative', K)
        colors = bmap.mpl_colors
        width=1
        xcoord = np.arange(len(data[a]['ids'])*width, step=width)
        bars=[]

        xcoord = np.arange(len(data[a]['ids'])*width, step=width)
        for i, (c, ki) in enumerate(zip(colors,data[a]['assn'])):
            if i ==0:
                bars.append( plt.bar(left=xcoord, height=ki, width=width, color=c, edgecolor='none') )
            else:
                bottom = np.sum(data[a]['assn'][:i,:], axis=0)
                bars.append( plt.bar(left=xcoord, height=ki, width=width, color=c, bottom=bottom, edgecolor='none') )

        plt.title('Structure Plot, K='+str(K))
        plt.grid(False)
        plt.xlim([0,xcoord[-1]+width])
        plt.ylim([0,1])
        plt.xticks(xcoord+width/2., data[a]['ids'], rotation=90, fontsize=6)
        plt.tick_params(axis='x', which='both', bottom='on', top='off', direction='out')
        plt.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')

        plt.show()


data = {}
for struct_file in iglob(getcwd()+'/*_f'):
        out = open(struct_file[:-2]+'_indiv_q.txt', 'w')
        started = False
        ids = []
        assignments = []
        sfn = open(struct_file)
        for l in sfn:
            if l == 'Inferred ancestry of individuals:\n':
                started = True
                sfn.next()
                continue
            if started and l == '\n':
                started = False
                out.close()
            if started:
                tmp = l.split()
                ids.append(tmp[1])
                assignments.append(tmp[5:])
                out.write('\t'.join(tmp[0:4])+'\t'+'\t'.join(tmp[5:])+'\n')
        data[struct_file[:-2]] = {'ids':ids, 'assn':np.array(assignments, dtype='float').T}

struct_to_fig(data)

