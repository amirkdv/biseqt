import os
from itertools import combinations
import logging
from biseqt.util import ProgressIndicator
from biseqt.blot import HomologyFinder
from biseqt.sequence import Alphabet
from biseqt.stochastics import rand_seq, MutationProcess, rand_read
from biseqt.kmers import KmerCache
from util import plot_classifier, log, with_dumpfile
from util import DATA_DIR, DUMP_DIR


# M is the mutation process to insert noise to each read, note that when
# comparing two reads later they have both suffered mutations at the rate
# dictated by M and hence are "twice" furhter apart as they are from the
# original genome.
def create_simulated_reads(M, path, gap=None, subst=None,
                           genome_size=10000, read_len=2000, num_reads=40):
    log('simulating %d sequencing reads of %d bp for genome of %d bp.' %
        (num_reads, read_len, genome_size))
    A = Alphabet('ACGT')
    M = MutationProcess(A, subst_probs=subst, ge_prob=gap, go_prob=gap)
    genome = rand_seq(A, genome_size)
    indic = ProgressIndicator(num_total=num_reads)
    indic.start()
    reads, mappings = [], []
    for read, start in rand_read(genome, len_mean=read_len, num=num_reads):
        indic.progress()
        reads.append(M.mutate(read)[0])
        mappings.append((start, start + len(read)))
    indic.finish()

    path = os.path.join(DATA_DIR, path)
    with open(path, 'w') as f:
        for m, read in zip(mappings, reads):
            f.write('> + %d:%d\n%s\n' % (m[0], m[1], str(read)))


def load_mapped_reads(path, max_num=100):
    A = Alphabet('ACGT')
    reads, mappings = [], []
    log('loading sequencing reads from %s' % path)
    count = 0
    reads, mappings = [], []
    with open(path) as f:
        while True:
            line = f.readline()
            if not line:
                break
            rc, coordinates = line[2:].split()
            assert rc in '+-'
            m0, m1 = coordinates.split(':')
            mappings.append((rc, int(m0), int(m1)))
            reads.append(A.parse(f.readline().strip()))
            count += 1
            if count == max_num:
                break
    return reads, mappings


@with_dumpfile
def sim_overlap_discovery(reads_path, wordlen=None, db_path=None, gap=None,
                          p_match=None, **kw):
    sim_data = {
        'scores': {
            's0_pos': [],
            's0_neg': [],
            's1_pos': [],
            's1_neg': [],
        },
        'gap': gap,
        'p_match': p_match,
        'wordlen': wordlen,
    }
    A = Alphabet('ACGT')
    db_path = os.path.join(DUMP_DIR, db_path)

    KC_kw = {'alphabet': A, 'wordlen': wordlen, 'path': db_path,
             'log_level': logging.WARN}
    kmer_cache = KmerCache(**KC_kw)
    HF_kw = {'g_max': gap, 'sensitivity': .9, 'kmer_cache': kmer_cache}
    HF_kw.update(**KC_kw)

    reads, mappings = load_mapped_reads(reads_path)

    def _overlaps(map0, map1):
        rc0, from0, to0 = map0
        rc1, from1, to1 = map1
        if rc0 != rc1:
            return False
        overlap_len = min(to0, to1) - max(from0, from1)
        if overlap_len > 10:
            return True
        else:
            return False

    num_reads = len(reads)
    assert len(reads) == len(mappings)
    num_total = (num_reads * (num_reads - 1)) / 2
    log('finding all overlapping pairs of reads')
    indic = ProgressIndicator(num_total=num_total, percentage=False)
    indic.start()
    for i, j in combinations(range(num_reads), 2):
        indic.progress()
        HF = HomologyFinder(reads[i], reads[j], **HF_kw)
        (seg0, s0), (seg1, s1) = HF.highest_scoring_overlap_band(p_min=p_match)
        if _overlaps(mappings[i], mappings[j]):
            sim_data['scores']['s0_pos'].append(s0)
            sim_data['scores']['s1_pos'].append(s1)
        else:
            sim_data['scores']['s0_neg'].append(s0)
            sim_data['scores']['s1_neg'].append(s1)

    indic.finish()
    return sim_data


def exp_overlap_detection():
    # NOTE when mutation probabilities for HF are lower than that of actual
    # mutaions we don't see any of the overlaps! (which is good). Also note
    # that since we are mutation both reads the HF and mutation probabilities
    # are significantly different.

    wordlen = 8
    hf_gap = .15
    hf_subst = .15
    hf_match = (1 - hf_gap) * (1 - hf_subst)

    reads_path = 'data/reads_sim.fa'
    db_path = 'overlaps.db'
    dumpfile = 'overlap_scores.txt'

    sim_data = sim_overlap_discovery(
        reads_path, wordlen=wordlen, db_path=db_path,
        gap=hf_gap, p_match=hf_match,
        dumpfile=dumpfile, ignore_existing=False)

    scores = sim_data['scores']
    labels = ['overlapping', 'non-overlapping']
    plot_classifier('overlap[H0].png', scores['s0_pos'], scores['s0_neg'],
                    labels=labels,
                    title='Overlap discovery with H0 (sim. data)')
    plot_classifier('overlap[H1].png', scores['s1_pos'], scores['s1_neg'],
                    labels=labels,
                    title='Overlap discovery with H1 (sim. data)')

    # ===============
    # Biological Data
    # ===============
    reads_path = 'data/leishmania/reads_mapped.fa'
    db_path = 'overlaps[bio].db'
    dumpfile = 'overlap_scores[bio].txt'

    sim_data = sim_overlap_discovery(
        reads_path, wordlen=wordlen, db_path=db_path,
        gap=hf_gap, p_match=hf_match,
        dumpfile=dumpfile, ignore_existing=False)

    scores = sim_data['scores']
    labels = ['overlapping', 'non-overlapping']
    plot_classifier('overlap[H0][bio].png', scores['s0_pos'], scores['s0_neg'],
                    labels=labels,
                    title='Overlap discovery with H0 (bio. data)')
    plot_classifier('overlap[H1][bio].png', scores['s1_pos'], scores['s1_neg'],
                    labels=labels,
                    title='Overlap discovery with H1 (bio. data)')


if __name__ == '__main__':
    # mut_gap = .08
    # mut_subst = .05
    # create_simulated_reads(M, 'reads_sim.fa', gap=mut_gap, subst=mut_subst)

    exp_overlap_detection()
