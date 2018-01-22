#!/usr/bin/env python
import pickle
from math import sqrt
from itertools import combinations

from biseqt.sequence import Alphabet
from biseqt.stochastics import MutationProcess
from biseqt.plots import figure, plot_roc, plot_cdf, save, plot_ppv, plot_npv
from biseqt.mapping import Read, ReadMapper
from biseqt.util import ProgressIndicator
from biseqt.stochastics import expected_alignment_length

from biseqt.seeds import SeedIndex

gap_prob = .15
wordlen = 10

A = Alphabet('ACGT')
subst_prob = .1
M = MutationProcess(A, subst_probs=subst_prob, go_prob=gap_prob, ge_prob=gap_prob)
subst_scores, (go_score, ge_score) = M.log_odds_scores()

#for mapping to reference
#A = Alphabet('ACGTN')
#subst_scores = [row + [0] for row in subst_scores] + [[0, 0, 0, 0, 0]]

min_band_score = -10
sensitivity = 1-1e-4
bin_len = 30
bands_kw = {
    'sensitivity': 1-1e-4,
    'gap_prob': gap_prob,
    #'max_kmer_score': 10,
    'max_kmer_score': 100,
}
aligner_kw = {
    #'max_new_mins': 10,
    #'min_score': 100,
    'subst_scores': subst_scores,
    'go_score': go_score,
    'ge_score': ge_score,
}

sam_path = 'leishmania/blasr.mappings.sam'
reads_fa = 'leishmania/reads.fa'
#refs_fa = 'leishmania/reference.upper.fa'
refs_fa = None

from multiprocess.managers import BaseManager
class LogManager(BaseManager):
    pass
LogManager.register('ProgressIndicator', ProgressIndicator)

from multiprocess import Manager
from multiprocess import Pool

# =========== Alignment classifier ====================
def aln_classifier_worker(args):
    mapper, indic, r0, r1, kw = args
    indic.progress()

    ours, target, aln= mapper.map_read(r0, r1, **kw)
    if ours is None or target is None:
        return {'ids': (r0.record.id, r1.id), 'aln': None}
    return {'ids': (ours.id, target.id), 'aln': aln}

def aln_classifier(mapper, overlaps):
    mapper.log('aligning all pairs')
    kw = {
        'gap_prob': gap_prob,
        'sensitivity': sensitivity,
        'bin_len': bin_len,
        #'min_band_score': min_band_score,
    }
    kw.update(aligner_kw)
    reads = mapper.load_reads()

    #m = LogManager()
    #m.start()
    #indic = m.ProgressIndicator(num_total=(len(reads) * (len(reads) - 1)), percentage=True)
    #indic.start()

    #def _pairs():
        #for r0, r1 in combinations(reads, 2):
            #yield mapper, indic, r0, r1.record, kw
            #yield mapper, indic, r0, r1.rc_record, kw

    ##res = []
    ##for args in _pairs():
        ##res.append(aln_classifier_worker(list(args)))

    #pool = Pool()
    #res = pool.map_async(aln_classifier_worker, _pairs())
    #pool.close()
    #pool.join()
    #res = res.get()

    #indic.finish()
    #mapper.log('Done mapping all reads')

    #res = {r['ids']: r['aln'] for r in res}
    #with open('alns.txt', 'w') as f:
        #pickle.dump(res, f)
        #mapper.log('dumped scores to %s' % f.name)

    with open('alns.txt') as f:
        res = pickle.load(f)

    pos_scores = []
    neg_scores = []
    for ids, aln in res.items():
        if aln is not None:
            pct_id = sum(float(op == 'M') for op in aln.transcript) / len(aln.transcript)
        else:
            pct_id = 0

        if aln is not None and ids in overlaps:
            pos_scores.append(pct_id)
        else:
            neg_scores.append(pct_id)

    plot_classifier(pos_scores, neg_scores, 'alignment', '', mapper.log, classifier='>')

#def plot_aln_classifier(mapper, overlaps):
    #max_id = 79
    #def _record_scores(f, pos_scores, neg_scores):
        #for line in f.readlines():
            #line = line.strip()
            #if not line:
                #continue
            #line = line.split(': ')
            #ids = eval(line[0])
            #if max(ids) > max_id:
                #continue
            #score, s, t, opseq = line[1].split(', ')
            #pct_id = sum(float(op == 'M') for op in opseq)/ len(opseq)
            #if ids in overlaps:
                #pos_scores.append(pct_id)
            #else:
                #neg_scores.append(pct_id)

    #pos_scores = []
    #neg_scores = []
    #with open('pos_aln.txt') as f:
        #_record_scores(f, pos_scores, neg_scores)
    #with open('neg_aln.txt') as f:
        #_record_scores(f, pos_scores, neg_scores)

    #print pos_scores
    #print '----------'
    #print neg_scores
    #plot_classifier(pos_scores, neg_scores, 'aln', '', mapper.log)

def plot_classifier(pos_scores, neg_scores, name, prefix, log, classifier='>'):
    # ---------- ROC -----------
    fig = figure(figsize=(12, 7))
    ax = fig.add_subplot(1, 2, 1)
    plot_cdf(ax, pos_scores, color='g', label='overlapping')
    plot_cdf(ax, neg_scores, color='r', label='non-overlapping')
    ax.set_title('Score CDF', fontsize=14)
    ax.legend(fontsize=12, loc='lower right')
    ax.tick_params(labelsize=14)

    ax = fig.add_subplot(1, 2, 2)
    ax.set_title('ROC', fontsize=14)
    plot_roc(ax, pos_scores, neg_scores, classifier=classifier)
    ax.tick_params(labelsize=14)

    path = '%s_%s_score_roc.png' % (prefix, name)
    save(fig, path, produce_legend=False, dpi=500)
    log('Saved ROC for %s classifier to %s' % (name, path))

    # -------- PPV -----------
    fig = figure()
    ax = fig.add_subplot(3, 1, 1)
    ax.set_xlabel('Band Score')
    plot_cdf(ax, pos_scores, color='g')
    plot_cdf(ax, neg_scores, color='r')
    plot_ppv(fig.add_subplot(3, 1, 2, sharex=ax), pos_scores, neg_scores, classifier=classifier)
    plot_npv(fig.add_subplot(3, 1, 3, sharex=ax), pos_scores, neg_scores, classifier=classifier)

    path = '%s_%s_predictive_value.png' % (prefix, name)
    save(fig, path)
    log('Saved predictive values for %s classifier to %s' % (name, path))

# ============= Both classifiers ============
def two_classifier_worker(args):
    mapper, indic, r0, r1, kw = args
    indic.progress()

    len0, len1 = r0.record.attrs['length'], r1.attrs['length']
    seed_index, ((id0, id1), band, score) = mapper.map_to_band(r0, r1,
        min_band_score=kw.get('min_band_score', None), gap_prob=kw['gap_prob'],
        sensitivity=kw['sensitivity'])
    assert band != (None, None)

    min_num = seed_index.min_num_seeds_in_band(id0, id1, len0, len1, bin_len=kw['bin_len'], diag_range=band)
    #print id0, id1, round(score, 2), min_num
    #print {'ids': (id0, id1), 'lens': (len0, len1), 'scores': (round(score, 2), min_num) }
    return {'ids': (id0, id1), 'lens': (len0, len1), 'scores': (round(score, 2), min_num) }

def two_classifier(mapper, prefix, overlaps):
    mapper.log('separating + and - sample sets by two classifiers')

    kw = {
        'gap_prob': gap_prob,
        'sensitivity': sensitivity,
        'bin_len': bin_len,
        #'min_band_score': min_band_score,
    }

    reads = mapper.load_reads()

    m = LogManager()
    m.start()
    indic = m.ProgressIndicator(num_total=(len(reads) * (len(reads) - 1)), percentage=True)
    indic.start()

    #for r0, r1 in combinations(reads, 2):
        #two_classifier_worker([mapper, indic, r0, r1.record, kw])
        #two_classifier_worker([mapper, indic, r0, r1.rc_record, kw])
    #raise

    pool = Pool()

    def _pairs():
        for r0, r1 in combinations(reads, 2):
            yield mapper, indic, r0, r1.record, kw
            yield mapper, indic, r0, r1.rc_record, kw

    res = pool.map_async(two_classifier_worker, _pairs())
    pool.close()
    pool.join()


    res = res.get()
    indic.finish()
    mapper.log('Done mapping all reads')

    res = {r['ids']: {
        'lens': r['lens'],
        'overlap': overlaps.get(r['ids'], (False, False))[0],
        'quality': overlaps.get(r['ids'], (False, False))[1],
        'scores': r['scores'],
    } for r in res}
    with open('%s_%d_double_scores.txt' % (prefix, bin_len), 'w') as f:
        pickle.dump(res, f)
        mapper.log('dumped scores to %s' % f.name)
    return

    cs, xs, ys = [], [], []
    first_pos_scores, second_pos_scores = [], []
    first_neg_scores, second_neg_scores = [], []
    for result_set in res:
        (id0, id1), scores = result_set['ids'], result_set['result']
        xs.append(scores[0])
        ys.append(scores[1])
        if (id0, id1) in overlaps:
            first_pos_scores.append(scores[0])
            second_pos_scores.append(scores[1])
            cs.append('g')
        else:
            first_neg_scores.append(scores[0])
            second_neg_scores.append(scores[1])
            cs.append('r')

    max_first = max(first_pos_scores)
    first_pos_scores = [s if s < float('inf') else max_first for s in first_pos_scores]

    plot_classifier(first_pos_scores, first_neg_scores, 'band', prefix, mapper.log, classifier='>')
    plot_classifier(second_pos_scores, second_neg_scores, 'alt-neg', prefix, mapper.log, classifier='<')

    fig = figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(xs, ys, color=cs, s=1, alpha=.7)
    ax.set_xlabel('Score in first classifier')
    ax.set_ylabel('Score in second classifier')

    path = '%s_double_scores.png' % prefix
    save(fig, path)
    mapper.log('Saved scores for to classifiers to ' + path)

if __name__ == '__main__':
    M = ReadMapper(A, wordlen, 'mem3_reads_test.db')

    from matplotlib import pyplot as plt
    ax = plt.axes()
    seed_index = SeedIndex(M.kmer_index)
    ids = (71, 94)
    seed_index.index_seeds(ids)
    seed_index.plot_seeds(ax, *ids, color='r')
    ax.set_ylim(12e3, 0)
    ax.set_xlim(0, 12e3)
    ax.set_xlabel('Position in read #%d' % ids[0])
    ax.set_ylabel('Position in read #%d' % ids[1])
    ax.grid(True)
    plt.savefig('prez-plots/fp.png', dpi=300)
    raise

    #M.initialize(reads_fa, refs_fa, num_reads=500)

    known_overlaps = list(M.overlaps_from_sam_mappings(sam_path))
    known_overlaps_by_id = {(_r1.id, _r2.id): (_len, _q) for _r1, _r2, _len, _q in known_overlaps}
    #known_overlaps_by_id = []

    #two_classifier(M, 'mem3', known_overlaps_by_id)
    #plot_aln_classifier(M, known_overlaps_by_id)
    aln_classifier(M, known_overlaps_by_id)
    raise

    #map_to_ref('ref')
    #fig = figure()
    #scores = [score for _, _, score in M.kmer_index.kmers()]
    #max_score = max(score for score in scores if score < float('+inf'))
    #scores = [scores if score < float(+'inf') else max_score for score in scores]
    #plot_cdf(fig.add_subplot(111), scores)
    #M.seed_index.plot_seeds(fig.add_subplot(111), 1, 21)
    #save(fig, 'debug.png')
