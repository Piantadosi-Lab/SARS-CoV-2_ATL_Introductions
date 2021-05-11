#!/usr/bin/env python3
import argparse
import baltic as bt
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import glob
import ast
from treetime.utils import numeric_date, datetime_from_numeric
import pandas as pd
import datetime
from collections.abc import Iterable
import json
import math
from utils import plot_style
import string
import numpy as np


def process_results(file, results):
    file_name = file.split('/')[-1].split('.')[0]
    dat = pd.read_csv(file, index_col = 0)
    dat = pd.melt(dat, id_vars=['destination region', 'source region', 'descendant tips', 'source node', '# descendants'], 
            value_vars=['midroot time', 'lower midroot time', 'upper midroot time'])
    dat['file'] = file_name
    dat = dat.sort_values('value')
    dat = dat.assign(cumsum=dat.groupby(['destination region', 'variable']).cumcount()+1)
    results = results.append(dat)
    return(results.reset_index(drop=True))



def summarize_results(results):
    qs = (0.025, 0.5, 0.975)
    regions = list(set(results['destination region']))
    variables = list(set(results['variable']))
    n_import_data = pd.DataFrame()
    observations = list(set(results['file']))
    # this is hacky, i'm bad at pandas groupby
    # we don't care about introductions leading to general wisconsin tips
    for region in ['MHD', 'UWHC']:  
        these_results = results[(results['destination region'] == region)]
        time_range = [min(these_results['value']), max(these_results['value'])]
        times = np.arange(np.floor(time_range[0]*100)/100, np.ceil(time_range[1]*100)/100, 1/1000)
        for variable in variables:
            these_results = results[(results['destination region'] == region) & 
                                    (results['variable'] == variable)]
            chunked_data = pd.DataFrame(index=times)
            chunked_data['region'] = region
            chunked_data['variable'] = variable
            chunked_data['importations'] = [[len(these_results[(these_results['value'] < time) & 
                                (these_results['file'] == observation)]) for observation in observations] for time in times]
            chunked_data['importations_summary'] = chunked_data['importations'].apply(hpd, qs=(0.025, 0.5, 0.975))
            chunked_data[['low', 'mid', 'high']] = \
                pd.DataFrame(chunked_data['importations_summary'].tolist(), index=chunked_data.index)
            n_import_data = n_import_data.append(chunked_data)
    return(n_import_data)


def plot_timeseries(ax, results, plot_config):
    [ax.spines[loc].set_visible(False) for loc in ['top', 'right']]
    ax.set_ylabel('Introductions')
    for idx, res in results.groupby(['file', 'destination region', 'variable']):
        res = res.sort_values('value')
        res_filtered = res[['value', 'cumsum']]
        res_filtered = pd.concat([res_filtered, 
            pd.DataFrame([[float(plot_config['max_date']), res_filtered['cumsum'].iloc[-1]]],
                columns=['value', 'cumsum'])])
        _ = ax.step(res_filtered['value'],
                 res_filtered['cumsum'], alpha=0.2, color=plot_config['colors'][idx[1]])
    y_lims= ax.get_ylim()
    return(ax, y_lims)



# //// HELPER FUNCTIONS FOR TREE PLOTTING ////
def branch_w_func(k, n_objects, branch_widths):
    if k.branchType == 'node':
        n_children = len(k.leaves)
        if n_children > n_objects*0.01: 
            return(branch_widths[0])
        else:
            return(branch_widths[1])
    else:
        return(branch_widths[2])


def branch_c_func(k, colors_dict, tree_states_dict):
    if k.branchType == 'node':
        if 'label' in k.traits.keys():
            name = k.traits['label']
        else:
            return(colors_dict['base'])
    else:
        name = k.name
    if tree_states_dict[name] in colors_dict.keys():
        return(colors_dict[tree_states_dict[name]])
    else:
        return(colors_dict['base'])


def tip_c_func(k, colors_dict, tree_states_dict):
    if tree_states_dict[k.name] in colors_dict.keys():
        return(colors_dict[tree_states_dict[k.name]])
    else:
        return('#ffffff')


def tip_s_func(k, tip_size, colors_dict, tree_states_dict):
    if tree_states_dict[k.name] in colors_dict.keys():
        # todo base this off axes size
        return(tip_size)
    else:
        return(0)


def prioritize_node(tree, node_id, prioritize=True):
    # hacky function to prioritze (y-axis position) of a specific node
    # get node in tree
    interest_node = \
        [i for i in tree.Objects if 
            i.branchType=='node' and i.traits['label'] == node_id][0]
    print(interest_node)
    # get parent of interest node
    interest_parent = interest_node.parent
    # get all children of parent node
    curr_children = interest_parent.children
    # remove interest node from list
    curr_children.remove(interest_node)
    if prioritize == True:
        # reinserts interest node at top of list
        curr_children.insert(0, interest_node)
        # replace children
        interest_parent.children = curr_children
    if prioritize == False:
        # reinserts interest node at top of list
        curr_children.insert(len(curr_children), interest_node)
        # replace children
        interest_parent.children = curr_children
    # redraw tree
    tree.drawTree()
    return(tree)


def import_tree(tree_file, max_date, node_name_len, prioritize_nodes=None, deprioritize_nodes=None):
    myTree = bt.loadNewick(tree_file)    
    myTree.setAbsoluteTime(max_date)
    # hack to remove bootstrap support values from node names
    for k in myTree.Objects:
        if k.branchType == 'node':
            if 'label' in k.traits.keys():
                k.traits['label'] = k.traits['label'][0:node_name_len]
    if prioritize_nodes:
        for node in prioritize_nodes:
            myTree = prioritize_node(myTree, node)
    if deprioritize_nodes:
        for node in deprioritize_nodes:
            myTree = prioritize_node(myTree, node, prioritize=False)
    #myTree.drawTree()
    return(myTree)


def plot_tree(ax, tree_file=None, tree_states_file=None, travel_file=None, plot_config=None, branch_widths=[1.5, 0.1, 0.05]):    
    tree_states = pd.read_csv(tree_states_file, sep=',')
    tree_states_dict = \
        {i['name']: i['state'] for idx, i in tree_states.iterrows()}
    for item in ['prioritize_nodes', 'deprioritize_nodes']:
        if item not in plot_config.keys():
            plot_config[item] = None
    myTree = \
        import_tree(tree_file, plot_config['max_date'], 12, 
            prioritize_nodes=plot_config['prioritize_nodes'], 
            deprioritize_nodes=plot_config['deprioritize_nodes'])
    xMin = min([i.absoluteTime for i in myTree.Objects])
    xMax = max([i.absoluteTime for i in myTree.Objects])
    xSpan = xMax - xMin
    L=len(list(filter(lambda k: k.branchType=='leaf',myTree.Objects)))
    tip_size = 50+myTree.treeHeight*35
    myTree.plotTree(ax, 
        x_attr=lambda k: k.absoluteTime,
        width=lambda k: branch_w_func(k, len(myTree.Objects), branch_widths), 
        colour=lambda k: branch_c_func(k, plot_config['colors'], tree_states_dict), 
        alpha=1.0, 
        connection_type='elbow')
    ax.set_xlim((xMin-xSpan*0.1, xMax+xSpan*0.1))
    myTree.plotPoints(ax, 
        x_attr=lambda k: k.absoluteTime, 
        size=lambda k: tip_s_func(k, tip_size, plot_config['colors'], tree_states_dict), 
        colour=lambda k: tip_c_func(k, plot_config['colors'], tree_states_dict),
        outline_colour=plot_config['colors']['base'])
    if travel_file != None:
        travel_tips = set(pd.read_csv(travel_file, header=None)[0])
        myTree.plotPoints(ax, 
        x_attr=lambda k: k.absoluteTime, 
        size=lambda k: tip_size if k.name.split('|')[1] in travel_tips else 0,
        colour=plot_config['colors']['traveler'],
        outline_colour=plot_config['colors']['base'], zorder=5)
    ax.set_ylim(-myTree.ySpan*0.01, myTree.ySpan*1.01)
    ax.set_yticks([])
    [ax.spines[loc].set_visible(False) for loc in ['top', 'right', 'left']]
    if 'labels' in plot_config.keys():
        for label_idx, (label_key, label) in enumerate(plot_config['labels'].items()):
            ax_height = (ax.get_position().y1 - ax.get_position().y0) * ax.figure.bbox.height
            font_size = 10+int(ax.figure.bbox.height*0.01)*2.5
            text = ax.text(-0.05, 0.15-(1.2*font_size/ax_height)*label_idx, 
                label, color=plot_config['colors'][label_key],
                        size=font_size, transform=ax.transAxes, 
                        path_effects=plot_config['effects'])
    return(ax, (xMin-xSpan*0.1, xMax+xSpan*0.1), myTree)



def add_heatmap(ax, tree=None, metadata_file=None, plot_config=None, max_val=1.0, width=0.01):
    metadata = pd.read_csv(metadata_file, sep='\t', header=None)
    clade_dict = \
        {plot_config['tree_name_dict'][i[1]]: i[8] for 
            idx, i in metadata.iterrows() if 
                i[1] in plot_config['tree_name_dict'].keys()}
    for k in tree.Objects:
        if k.branchType == 'leaf':
            # a small number of sequences do not have clades assigned
            if clade_dict[k.name] in plot_config['colors'].keys():
                ax.add_patch(Rectangle((max_val,k.y-0.5),width,1,
                    facecolor=plot_config['colors'][clade_dict[k.name]], 
                    edgecolor='none'))        
    # need to label clades
    # todo better (defaultdict)
    y_val_dict = \
        {k.name: k.y for k in tree.Objects if k.branchType == 'leaf'}
    all_clades = set([i for i in clade_dict.values() if i==i])
    for clade in all_clades:
        # get all seqs belonging to this clade
        clade_seqs = \
            [k for k in y_val_dict.keys() if clade_dict[k] == clade]
        # get y-values of all seqs
        clade_y = [y_val_dict[i] for i in clade_seqs]
        mean_clade_y = sum(clade_y)/len(clade_y)
        ax.text(max_val+width*1.10, mean_clade_y, clade, 
            horizontalalignment='left', verticalalignment='center', 
            size=6+int(ax.figure.bbox.height*0.01)*1.25,
            color=plot_config['colors'][clade], path_effects=plot_config['effects'])
    return(ax, max_val+width)


def plot_clade_distribution(ax, metadata_file, plot_config):
    metadata = pd.read_csv(metadata_file, sep='\t', header=None)
    metadata = metadata[metadata[6].isin(plot_config['colors'].keys())]
    metadata['week_end_date'] = pd.to_datetime(metadata[2]).apply(lambda k: k+datetime.timedelta(days= 6 - k.weekday()))
    n_per_week = metadata.groupby([8, 'week_end_date']).size().reset_index()
    n_per_week = n_per_week.pivot(index='week_end_date', columns=8, values=0).fillna(0)
    n_per_week.index = n_per_week.index.to_series().apply(numeric_date)
    cm = LinearSegmentedColormap.from_list('clades', [plot_config['colors'][i] for i in n_per_week.columns], N=n_per_week.shape[1])
    #(n_per_week.div(n_per_week.sum(axis=1), axis=0)*100).plot.area(ax=ax, colormap=cm, alpha=1.0, zorder=2)
    n_per_week.plot.area(ax=ax, colormap=cm, alpha=1.0, zorder=2)
    n_per_week['dt'] = n_per_week.index.to_series().apply(lambda k: datetime_from_numeric(k))
    n_per_week.to_csv('.'.join(metadata_file.split('.')[0:-1])+'_n_clade_per_week.csv')
    ax.set_ylabel('Sequences/Wk', size=12)
    ax.get_legend().remove()
    ax.set_xlabel(None)
    return(ax)





def run():
    parser = argparse.ArgumentParser()
    # input files
    parser.add_argument('--tree', 
                        help='newick tree file to plot',
                        default=None)
    parser.add_argument('--treeNameSep', 
                        help='delimiter to seperate tree names on to match with metadata',
                        default='|')
    parser.add_argument('--treeNameField', 
                        help='which field in split tree names to match to metadata',
                        default=1,
                        type=int)
    parser.add_argument('--treeBootstrapDir', 
                        help='directory with bootstrap replicate trees',
                        default=None)
    parser.add_argument('--treeStates',
        help='file with importation data from ML tree')
    parser.add_argument('--metadata',
        help='file with metadata')
    parser.add_argument('--travelData',
        help='file with metadata')
    parser.add_argument('--config',
        help='json file with config information (labels, colors, etc.)')
    args = parser.parse_args()
        #args.tree = 'data/weighted_downsampling/ga_focused_aligned_masked_weighted.treefile_tres/0/0_refined_time.newick'
        #args.treeStates = 'data/weighted_downsampling/ga_focused_aligned_masked_weighted.treefile_tres/0/0_refined_node_states.csv'
        #args.treeBootstrapDir = 'data/weighted_downsampling/ga_focused_aligned_masked_weighted.ufboot_tres/*/*_importations.csv'
        #args.config = 'config/annotated_tree.json'
        #args.metadata = 'data/weighted_downsampling/ga_focused_aligned_masked_weighted_all_included_seqs.tsv'
        # //// CONFIGURES SETTINGS FOR PLOT ////
    plot_style()
    plot_config = json.load(open(args.config, 'r'))
    plot_config['effects'] = eval(plot_config['effects'])
    if 'prioritize_nodes' in plot_config.keys():
        plot_config['prioritize_nodes'] = eval(plot_config['prioritize_nodes'])
    if 'deprioritize_nodes' in plot_config.keys():
        plot_config['deprioritize_nodes'] = eval(plot_config['deprioritize_nodes'])
    plot_config['max_date'] = float(plot_config['max_date'])
    # todo generate programmaticaly
    x_ticks = [datetime.datetime.strptime('2020-01-01', '%Y-%m-%d'),
        datetime.datetime.strptime('2020-01-15', '%Y-%m-%d'),
        datetime.datetime.strptime('2020-02-01', '%Y-%m-%d'),
        datetime.datetime.strptime('2020-02-15', '%Y-%m-%d'),
        datetime.datetime.strptime('2020-03-01', '%Y-%m-%d'),
        datetime.datetime.strptime('2020-03-15', '%Y-%m-%d'),
        datetime.datetime.strptime('2020-04-01', '%Y-%m-%d')]
    x_labels = \
        [i.strftime("%m")+'/'+i.strftime("%d") for i in x_ticks]
    x_ticks = [numeric_date(i) for i in x_ticks]
    # //// PROCESS BOOTSTRAP RESULTS ////
    results = pd.DataFrame()
    for file in glob.glob(args.treeBootstrapDir):
        results = process_results(file, results)
    print(results)
    midroot_introductions = \
        results[results['variable'] == 'midroot time'].groupby(
            ['file', 'destination region', 'variable'])['cumsum'].agg('max').reset_index()
    # //// READ IN AND PROCESS ML TREE TREE STATES////
    # Figure 1
    fig = plt.figure(figsize=(6.4,4.8*2.5))
    plt.rcParams['hatch.color'] = plot_config['colors']['base']
    #plt.rcParams['hatch.linewidth'] = 0.5
    gs=GridSpec(7,16, figure=fig)
    ax00 = fig.add_subplot(gs[0,:])
    ax0 = fig.add_subplot(gs[1:6,:])
    ax1 = fig.add_subplot(gs[6,:])
    ax2 = fig.add_subplot(gs[6,14:])
    ax00 = plot_clade_distribution(ax00, args.metadata, plot_config)
    [ax00.spines[loc].set_visible(False) for loc in ['top', 'right']]



    ax0, x_lims, myTree = plot_tree(ax0, tree_file=args.tree, travel_file=args.travelData, 
        tree_states_file=args.treeStates, 
        plot_config=plot_config, branch_widths=[2, 0.02, 0.01])
    tree_name_dict = \
       {i.name.split(args.treeNameSep)[args.treeNameField]: i.name for 
           i in myTree.Objects if i.branchType == 'leaf'}
    plot_config['tree_name_dict'] = tree_name_dict
    ax0, x_max = \
        add_heatmap(ax0, tree=myTree, metadata_file=args.metadata, plot_config=plot_config, max_val=x_lims[1])

    x_lims = [x_lims[0], x_lims[1] + 7/366]
    ax1, y_lims = plot_timeseries(ax1, results, plot_config)


    sns.kdeplot(midroot_introductions["cumsum"], ax=ax2, shade=True, color='#8FBCBB', alpha=0.85, vertical=True)
    [ax2.spines[loc].set_visible(False) for loc in ['top', 'right', 'left', 'bottom']]
    ax2.patch.set_alpha(0)
    ax2.set_yticks([])
    ax2.set_xticks([])
    ax2.set_ylabel(None)
    ax2.set_xlabel(None)
    ax2.set_ylim(y_lims)
    for ax in [ax00, ax0, ax1]: 
        ax.set_xlim(x_lims)
        ax.set_xticks(x_ticks)
        ax.grid(axis='x', ls='-', color='#d8dee9', zorder=0)


    ax00.set_xticklabels([])
    ax0.set_xticklabels([])
    ax1.set_xticklabels(x_labels)
    ax1.set_xlabel('Date (2020)')
    fig.tight_layout()

    y_poses = [1.0, 1.0, 1.08]
    for ax_idx, ax in enumerate([ax00, ax0, ax1]):
        ax.text(-0.18, y_poses[ax_idx], string.ascii_uppercase[ax_idx], transform=ax.transAxes, 
                size=20, weight='bold', va="top")

    fig.tight_layout()
    #fig.subplots_adjust(hspace=0.05)
    fig.savefig(plot_config['out_name']+'.pdf')
    plt.close()



if __name__ == "__main__":
    run()
