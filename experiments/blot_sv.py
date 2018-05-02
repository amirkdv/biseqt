#!/usr/bin/env python
import logging
from itertools import combinations, product
import sys
import numpy as np
import igraph
from time import time
from matplotlib import pyplot as plt

from biseqt.blot import WordBlotLocalRef
from biseqt.sequence import Alphabet
from biseqt.stochastics import rand_seq, rand_read, MutationProcess

from util import plot_with_sd, savefig, with_dumpfile, log


def gen_data_set(**kw):
    alphabet, sv_type = kw['alphabet'], kw['sv_type']
    ref_len, sv_len, sv_pos = kw['ref_len'], kw['sv_len'], kw['sv_pos']
    gap, subst = kw['gap'], kw['subst']
    coverage, margin = kw['coverage'], kw['margin']
    sequencing_kw = {
        'len_mean': sv_len * 2,
        'len_sd': sv_len * .1,
        'expected_coverage': coverage,
    }

    ref = rand_seq(alphabet, ref_len)
    if sv_type == 'insertion':
        indiv = ref[:sv_pos] + rand_seq(alphabet, sv_len) + ref[sv_pos:]
    elif sv_type == 'deletion':
        indiv = ref[:sv_pos] + ref[sv_pos + sv_len:]
    elif sv_type == 'duplication':
        sv_src = kw['sv_src'], kw['sv_src'] + sv_len
        indiv = ref[:sv_pos] + ref[sv_src[0]:sv_src[1]] + ref[sv_pos:]
    else:
        raise ValueError('unknown SV type %s' % sv_type)

    M = MutationProcess(alphabet, subst_probs=subst, ge_prob=gap, go_prob=gap)
    reads = [(M.mutate(read)[0], pos)
             for read, pos in rand_read(indiv, **sequencing_kw)]
    labels = [False] * len(reads)
    for idx, (read, start_pos) in enumerate(reads):
        overlap_with_sv = min(start_pos + len(read), sv_pos + sv_len) - \
                          max(start_pos, sv_pos)
        if overlap_with_sv < margin:
            continue

        start_overlap = start_pos < sv_pos
        end_overlap = start_pos + len(read) > sv_pos + sv_len
        if sv_type == 'insertion':
            labels[idx] = start_overlap and end_overlap
        elif sv_type == 'deletion':
            labels[idx] = start_overlap
        elif sv_type == 'duplication':
            labels[idx] = start_overlap or end_overlap
        else:
            raise ValueError('unknown SV type %s' % sv_type)

    reads = [read for read, pos in reads]
    return {'ref': ref, 'indiv': indiv, 'reads': reads, 'labels': labels}


def chain_paths_on_read(read, sims, **kw):
    margin = kw['margin']
    seg_graph = igraph.Graph()
    for sim in sims:
        seg_graph.add_vertex(rec=sim)
    seg_graph.add_vertex(name='start')
    seg_graph.add_vertex(name='end')
    for idx, sim in enumerate(sims):
        if sim['read'][0] < margin:
            seg_graph.add_edge('start', idx)
        if sim['read'][1] > len(read) - margin:
            seg_graph.add_edge(idx, 'end')

    for (idx0, sim0), (idx1, sim1) in combinations(enumerate(sims), 2):
        from0, to0 = sim0['read']
        from1, to1 = sim1['read']
        overlap_len = min(to0, to1) - max(from0, from1)
        if overlap_len > -margin:
            # HACK what to do with direction?
            if from0 <= from1:
                seg_graph.add_edge(idx0, idx1)
            else:
                seg_graph.add_edge(idx1, idx0)
    return seg_graph.get_all_shortest_paths('start', to='end', mode='out')


def chain_paths_on_ref(read, sims, **kw):
    margin = kw['margin']
    seg_graph = igraph.Graph()
    for sim in sims:
        seg_graph.add_vertex(rec=sim)
    seg_graph.add_vertex(name='start')
    seg_graph.add_vertex(name='end')
    min_start = min(sim['ref'][0] for sim in sims)
    max_end = max(sim['ref'][1] for sim in sims)
    for idx, sim in enumerate(sims):
        if sim['ref'][0] == min_start:
            seg_graph.add_edge('start', idx)
        if sim['ref'][1] == max_end:
            seg_graph.add_edge(idx, 'end')

    for (idx0, sim0), (idx1, sim1) in combinations(enumerate(sims), 2):
        from0, to0 = sim0['ref']
        from1, to1 = sim1['ref']
        overlap_len = min(to0, to1) - max(from0, from1)
        if overlap_len > -margin:
            # HACK what to do with direction?
            if from0 <= from1:
                seg_graph.add_edge(idx0, idx1)
            else:
                seg_graph.add_edge(idx1, idx0)
    return seg_graph.get_all_shortest_paths('start', to='end', mode='out')


# mode tells us whether the path is on read or on ref
def sv_coords(sims, path, mode='read'):
    assert len(path) == 4
    segs_in_order = sorted(path[1:-1], key=lambda i: sims[i]['ref'][0])
    coord0 = sims[segs_in_order[0]]['ref']
    coord1 = sims[segs_in_order[1]]['ref']
    if mode == 'read':
        return coord0[1], coord1[0]
    else:
        return (coord0[1] + coord1[0]) / 2


def call_svs(ref, reads, margin=50, K_min=200, p_min=.8, **WB_kw):
    """Calls structural variants based on single read mappings. For each read,
    the algorithm proceeds as follows:

    .. figure::
        https://www.dropbox.com/s/rjqobtkp60snte7/SV-schematic.png?raw=1
       :target:
        https://www.dropbox.com/s/rjqobtkp60snte7/SV-schematic.png?raw=1
       :alt: lightbox

       Schematic representation of the chain graphs used to call structural
       variants. Each read (vertical axis) is compared against the reference
       sequence (horizontal axis) and in each case two chain graphs are created
       on each of the read and reference axes where two similar segments are
       connected with an edge if their projections on the corresponding axis
       can be consistently chained.

    * all local similarities are found via Word-Blot
      (:func:`biseqt.blot.WordBlotLocalRef.similar_segments` with given
      :math:`K_{\min}, p_{\min}`).
    * Two chain graphs are built based on the projections of similarities on
      the read and on the reference genome, henceforth the *read graph* and the
      *reference graph*.
    * Reads containining SVs are identified using the following rules:
        * Normal reads are characterized by a shortest path in both the read
          and reference graphs with exactly two edges between the start and the
          end of projections.
        * deletions and duplications produce shortest paths with four edges in
          the read graph.
        * insertions produce shortest paths with four edges in the reference
          graph.
    """
    WB = WordBlotLocalRef(ref, **WB_kw)
    label_hats = [{'I': [False, []],  # insertion
                   'D': [False, []]}  # deletion / duplication
                  for _ in reads]

    for read_idx, read in enumerate(reads):
        sims = list(
            WB.similar_segments(read, K_min, p_min)
        )
        for sim_idx, sim in enumerate(sims):
            ref_pos, read_pos = WB.to_ij_coordinates_seg(sim['segment'])
            sims[sim_idx]['read'] = read_pos
            sims[sim_idx]['ref'] = ref_pos

        # for duplication or deletion
        d_paths = chain_paths_on_read(read, sims, margin=margin)
        if d_paths:
            d_path = d_paths[0]
            label_hats[read_idx]['D'][0] = len(d_path) > 3
            if len(d_path) == 4:
                label_hats[read_idx]['D'][1] = sv_coords(sims, d_path,
                                                         mode='read')

        i_paths = chain_paths_on_ref(read, sims, margin=margin)
        if i_paths:
            i_path = i_paths[0]
            label_hats[read_idx]['I'][0] = len(i_path) > 3
            if len(i_path) == 4:
                label_hats[read_idx]['I'][1] = sv_coords(sims, i_path,
                                                         mode='ref')

        sys.stderr.write('.')

    return label_hats


@with_dumpfile
def sim_structural_variants(**kw):
    sv_len, ps, wordlen = kw['sv_len'], kw['ps'], kw['wordlen']
    n_samples, coverage, margin = kw['n_samples'], kw['coverage'], kw['margin']
    sv_types = ['insertion', 'deletion', 'duplication']
    # maps sv types to sv call types
    sv_type_dict = {'insertion': 'I', 'duplication': 'D', 'deletion': 'D'}
    A = Alphabet('ACGT')
    WB_kw = {'g_max': .2, 'sensitivity': .9, 'alphabet': A, 'wordlen': wordlen,
             'log_level': logging.WARN}

    def _zero(): return np.zeros((len(ps), n_samples))

    sim_data = {
        'coverage': coverage,
        'n_samples': n_samples,
        'sv_len': sv_len,
        'ps': ps,
        'WB_kw': WB_kw,
        'coords': {sv_type: {p_idx: {idx: [] for idx in range(n_samples)}
                             for p_idx in range(len(ps))}
                   for sv_type in sv_types},
        'stats': {stat: {sv_type: _zero() for sv_type in sv_types}
                  for stat in ['tpr', 'fpr']},
        'times': _zero(),
    }

    for sv_type, (p_idx, p_match) in product(sv_types, enumerate(ps)):
        ref_len = sv_len * 5
        sv_src = sv_len
        sv_pos = sv_len * 3

        # distribute p_match evenly over gap and subst
        subst = gap = 1 - np.sqrt(p_match)

        log('SV: %s, p = %.2f (n=%d rounds)' %
            (sv_type, p_match, n_samples))
        for sample_idx in range(n_samples):
            dataset = gen_data_set(
                sv_type=sv_type,
                alphabet=A,
                ref_len=ref_len,
                sv_len=sv_len,
                sv_src=sv_src,
                sv_pos=sv_pos,
                gap=gap,
                subst=subst,
                coverage=coverage,
                margin=margin,
            )
            labels = dataset['labels']
            K_min = 200
            p_min = .6
            t_start = time()
            label_hats = call_svs(dataset['ref'], dataset['reads'],
                                  K_min=K_min, p_min=p_min, **WB_kw)
            sim_data['times'][p_idx, sample_idx] = time() - t_start

            tp, fp, tn, fn = 0, 0, 0, 0
            for label, label_hat in zip(labels, label_hats):
                key = sv_type_dict[sv_type]
                if label and label_hat[key][0]:
                    tp += 1
                elif label and not label_hat[key][0]:
                    fn += 1
                elif not label and label_hat[key][0]:
                    fp += 1
                elif not label and not label_hat[key][0]:
                    tn += 1

                if label_hat[key][0]:
                    sim_data['coords'][sv_type][p_idx][sample_idx].append(
                        label_hat[key][1])

            # almost imposible that tn + fp = 0, but tp + fn can be zero:
            if tp + fn == 0:
                tpr = 1.
            else:
                tpr = 1. * tp / (tp + fn)
            fpr = 1. * fp / (tn + fp)
            sim_data['stats']['tpr'][sv_type][p_idx, sample_idx] = tpr
            sim_data['stats']['fpr'][sv_type][p_idx, sample_idx] = fpr
            sys.stderr.write(' tpr = %.2f, fpr = %.2f' % (tpr, fpr))
            sys.stderr.write('\n')

    return sim_data


def plot_structural_variants(sim_data, suffix=''):
    ps = sim_data['ps']
    sv_types = ['duplication', 'deletion', 'insertion']

    fig = plt.figure(figsize=(11, 4))
    ax_stats = {}
    # ax_coord = {}
    for idx, sv_type in enumerate(sv_types):
        ax_stats[sv_type] = fig.add_subplot(1, 3, idx + 1)

    kw = {'markersize': 5, 'marker': 'o', 'alpha': .7}
    ls = {'tpr': '-', 'fpr': '--'}
    label = {'tpr': 'TPR (sensitivity)', 'fpr': 'FPR (1 - specificity)'}
    for sv_type in sv_types:
        for stat, values in sim_data['stats'].items():
            plot_with_sd(ax_stats[sv_type],
                         ps, sim_data['stats'][stat][sv_type][:, :], axis=1,
                         ls=ls[stat], label=label[stat],
                         **kw)
            ax_stats[sv_type].set_ylim(-.1, 1.3)
        ax_stats[sv_type].set_title(sv_type)
        ax_stats[sv_type].legend(loc='upper left')
        ax_stats[sv_type].set_xlabel('match probability')
        ax_stats[sv_type].set_xticks(ps)
        ax_stats[sv_type].set_xticklabels(ps)

    times = sim_data['times'].flatten()
    print 'time: %.2f s (%.2f)' % (times.mean(), np.sqrt(times.var()))

    savefig(fig, 'structural_variants%s.png' % suffix)


def exp_structural_variants():
    """Performance of Word-Blot in detecting structural variations (SV) in
    *simulations*:

    * Given a single read containing part of a SV and a reference sequence, one
      can deduce the presence of a SV from the topological arrangement of local
      similarities. The simple algorithm used here, :func:`call_svs`, simply
      pays attention to shortest paths in either of two directed graphs
      obtained by chaining projected similarities on either the reference or
      the read sequences.
    * For each of three copy number variation SVs, *insertion, deletion,
      duplication*, similar segments detected by Word-Blot are used to detect
      whether each read contains an SV. In each case the true positive rate and
      false positive rate are plotted as a function of similarity between
      individual and reference sequences.

    **Supported Claims**

    * The simple topological algorithm in :func:`call_svs` can accurately
      distinguish normal reads from those containing a SV.
    * The crux of our simple, single-read SV calling algorithm is that,
      roughly, SVs can be detected in a single read whenever there are more
      than one read/reference local similarities that, when chained, span the
      entire read (for duplication and deletion) or a whole region on reference
      (for insertion). This shows that Word-Blot accurately recovers local
      similarities at their maximal length without producing unnecessary
      fragmentation.
    * Word-Blot accurately recovers local similarities to their maximal extent
    * Since Word-Blot is a general purpose local similarity search and thus no
      assumption of collinearity it can accurately identify simple SVs (copy
      number variations) by simply looking at mapping of individual reads to a
      reference sequence.

    .. figure::
        https://www.dropbox.com/s/ase6z7wzlbc6dy0/
        structural_variants%5Bw%3D8%5D.png?raw=1
       :target:
        https://www.dropbox.com/s/ase6z7wzlbc6dy0/
        structural_variants%5Bw%3D8%5D.png?raw=1
       :alt: lightbox

       True positive rate (solid) and false positive rate (dashed) for
       discriminating normal reads from those containing evidence for SV,
       duplications (*left*), deletion (*middle*), insertion (*right*) using
       the simple algorithm based on Word-Blot local similarities. SV length is
       500nt, reference sequence length is 2500nt, sequencing coverage is 10,
       word length is 8, and each condition (match probability and SV type) is
       repeated n=5 times, average computation time for each read is 0.8
       seconds.
    """

    wordlen = 8  # kept in memory; don't go too high up
    sv_len = 500
    margin = 50
    ps = [round(.98 - .03 * i, 2) for i in range(6)]
    coverage = 10
    n_samples = 5

    suffix = '[w=%d]' % wordlen
    dumpfile = 'sv_calling%s.txt' % suffix

    sim_data = sim_structural_variants(
        sv_len=sv_len, ps=ps, wordlen=wordlen, n_samples=n_samples,
        margin=margin, coverage=coverage,
        dumpfile=dumpfile, ignore_existing=False,
    )
    plot_structural_variants(sim_data, suffix=suffix)


if __name__ == '__main__':
    exp_structural_variants()
