"""
Rountines for generating dotplots.
"""

import collections
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

import kanapy


def dotplot(seq_x, seq_y, config=dict(), title=None, anno_list=None):

    # Get k-mer util.
    if 'kutil' in config:
        kutil = config['kutil']

    if anno_list is None:
        anno_list = list()

    else:
        if 'k' in config:
            kutil = kanapy.util.kmer.KmerUtil(np.int32(config['k']))
        else:
            kutil = kanapy.util.kmer.KmerUtil(32)

    # Get sequence labels and start positions
    label_x = config.get('label_x', 'Sequence')
    start_x = config.get('start_x', 0)
    label_y = config.get('label_y', 'Reference')
    start_y = config.get('start_y', 0)

    linewidth = config.get('linewidth', 1)

    plot_width = config.get('plot_width', 7)
    plot_height = config.get('plot_height', 7)
    plot_dpi = config.get('plot_dpi', 300)
    invert_y = config.get('invert_y', False)

    plot_props = dict()  # Plot properties, such as x and y limits

    # Process annotation dict
    anno_list_background = list()
    anno_list_foreground = list()

    if anno_list is None:
        anno_dict = list()

    anno_element_index = 0

    for anno_element in anno_list:

        anno_element['index'] = anno_element_index

        if anno_element.get('background', True):
            anno_list_background.append(anno_element)
        else:
            anno_list_foreground.append(anno_element)

        anno_element_index += 1

    # Get Y k-mers and indexes in fwd and reverse orientation
    y_kmer_fwd = collections.defaultdict(set)
    y_kmer_rev = collections.defaultdict(set)

    for kmer, index in kanapy.util.kmer.stream(seq_y, kutil, True):
        y_kmer_fwd[kmer].add(index)
        y_kmer_rev[kutil.rev_complement(kmer)].add(index)

    # Get X k-mers
    x_mer_index = list(x_tuple for x_tuple in kanapy.util.kmer.stream(seq_x, kutil, True))

    # Get dotplot points and limits
    points_fwd = {(x + start_x, y + start_y) for kmer, x in x_mer_index for y in y_kmer_fwd[kmer]}
    points_rev = {(x + start_x, y + start_y) for kmer, x in x_mer_index for y in y_kmer_rev[kmer]}

    plot_props['min_x'] = min([
        min(points_fwd, key=lambda val: val[0])[0] if points_fwd else np.nan,
        min(points_rev, key=lambda val: val[0])[0] if points_rev else np.nan
    ])

    plot_props['max_x'] = max([
        max(points_fwd, key=lambda val: val[0])[0] if points_fwd else np.nan,
        max(points_rev, key=lambda val: val[0])[0] if points_rev else np.nan
    ])

    plot_props['min_y'] = min([
        min(points_fwd, key=lambda val: val[1])[1] if points_fwd else np.nan,
        min(points_rev, key=lambda val: val[1])[1] if points_rev else np.nan
    ])

    plot_props['max_y'] = max([
        max(points_fwd, key=lambda val: val[1])[1] if points_fwd else np.nan,
        max(points_rev, key=lambda val: val[1])[1] if points_rev else np.nan
    ])

    # Collapse points to lines - fwd
    lines_fwd = list()

    while points_fwd:
        x_start, y_start = min(points_fwd)
        points_fwd.discard((x_start, y_start))

        x_end = x_start + 1
        y_end = y_start + 1

        while (x_end, y_end) in points_fwd:
            points_fwd.discard((x_end, y_end))

            x_end += 1
            y_end += 1

        lines_fwd.append(((x_start, y_start), (x_end, y_end)))

    # Collapse points to lines - rev
    lines_rev = list()

    while points_rev:
        x_start, y_start = min(points_rev, key=lambda vals: (vals[0], -vals[1]))  # Min x, max y

        points_rev.discard((x_start, y_start))

        x_end = x_start + 1
        y_end = y_start - 1

        while (x_end, y_end) in points_rev:
            points_rev.discard((x_end, y_end))

            x_end += 1
            y_end -= 1

        lines_rev.append(((x_start, y_start), (x_end, y_end)))


    # Make plot
    fig = plt.figure(figsize=(plot_width, plot_height), dpi=plot_dpi)

    ax1 = fig.subplots(1, 1)

    # Background annotations
    for anno_element in anno_list_background:
       _plot_anno_element(ax1, anno_element, plot_props)

    # Dot plot points (as connected lines)
    for (x1, y1), (x2, y2) in lines_fwd:
        ax1.plot((x1, x2), (y1, y2), linewidth=linewidth, color='black')

    for (x1, y1), (x2, y2) in lines_rev:
        ax1.plot((x1, x2), (y1, y2), linewidth=linewidth, color='red')

    # Foreground annotations
    for anno_element in anno_list_foreground:
       _plot_anno_element(ax1, anno_element, plot_props)

    # Set x and y ranges
    ax1.set_xlim((start_x, start_x + len(seq_x)))
    ax1.set_ylim((start_y, start_y + len(seq_y)))

    # Invert y
    if invert_y:
        ax1.invert_yaxis()

    # Plot aestetics
    ax1.get_xaxis().set_major_formatter(
        mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ','))
    )

    ax1.get_yaxis().set_major_formatter(
        mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ','))
    )

    plt.xticks(rotation=45, fontsize=8, horizontalalignment='right')
    plt.yticks(rotation=45, fontsize=8, verticalalignment='top')

    fig.subplots_adjust(bottom=0.15)  # Test/tune

    if label_x is not None:
        ax1.set_xlabel(label_x)

    if label_y is not None:
        ax1.set_ylabel(label_y)

    if title is not None:
        fig.suptitle(title)

    mpl.rc('xtick', labelsize=5)

    # Return figure
    return fig

def _plot_anno_element(ax, anno_element, plot_props):

    anno_element_index = anno_element.get('index')

    if 'type' not in anno_element:
        raise RuntimeError(
            'Missing annotation type: "annotype" for annotation list element {}'.format(anno_element_index)
        )

    type = anno_element['type']

    if type == 'vshade':
        if 'x1' not in anno_element or 'x2' not in anno_element:
            raise RuntimeError(
                'Missing x1 and/or x2 for vshade annotation: element {}'.format(anno_element_index)
            )

        ax.axvspan(
            anno_element.get('x1'),
            anno_element.get('x2'),
            color=anno_element.get('color', 'black'),
            alpha=anno_element.get('alpha', 0.2)
        )

    elif type == 'hshade':
        if 'y1' not in anno_element or 'y2' not in anno_element:
            raise RuntimeError(
                'Missing y1 and/or y2 for hshade annotation: element {}'.format(anno_element_index)
            )

        x = np.array([plot_props['min_x'], plot_props['max_x']])
        y1 = np.repeat(anno_element.get('y1'), 2)
        y2 = np.repeat(anno_element.get('y2'), 2)

        ax.fill_between(
            x=x, y1=y1, y2=y2,
            color=anno_element.get('color', 'black'),
            alpha=anno_element.get('alpha', 0.2)
        )

    elif type == 'hlines' or type == 'hline':
        if 'y' not in anno_element:
            raise RuntimeError(
                'Missing y for hlines annotation: element {}'.format(anno_element_index)
            )

        ax.hlines(
            anno_element.get('y'),
            color=anno_element.get('color', 'black'),
            alpha=anno_element.get('alpha', 1.0),
            linestyles=anno_element.get('style', 'solid')
        )

    elif type == 'vlines' or type == 'vline':

        if 'x' not in anno_element:
            raise RuntimeError(
                'Missing x for vlines annotation: element {}'.format(anno_element_index)
            )

        if 'ymin' not in anno_element:
            raise RuntimeError(
                'Missing ymin for vlines annotation: element {}'.format(anno_element_index)
            )

        if 'ymax' not in anno_element:
            raise RuntimeError(
                'Missing ymax for vlines annotation: element {}'.format(anno_element_index)
            )

        ax.vlines(
            x=anno_element.get('x'),
            ymin=anno_element.get('ymin'),
            ymax=anno_element.get('ymax'),
            color=anno_element.get('color', 'black'),
            alpha=anno_element.get('alpha', 1.0),
            linestyles=anno_element.get('style', 'solid')
        )

    else:
        raise RuntimeError('Unknown annotation type: {}'.format(type))
