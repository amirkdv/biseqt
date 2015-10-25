#!/usr/bin/env python
from termcolor import colored
from . import scan

def print_alignment(S, T, transcript, f, width=120, margin=20, colors=True):
    """Deciphers a given transcript and writes pretiffied output to f.
    Transcripts are expected to have the following format:
        (<Si,Tj>),<score>:...

    :param S: "from" sequence (anything that casts to the correct str object).
    :param T: "to" sequence (like S).
    :param transcript: the transcript output of libalign, see
        `AlignProblem.solve`.
    :param f: file handle to write the output to; for standard output use
        `sys.stdout`.
    :param width(optional): terminal width (int) used for wrapping; default 120.
    :param margin(optional): length (int) of leading and trailing sequences in S
        and T before and after the alignment; default 20.
    :param colors(optional): whether or not (truthy) to use colors in output;
        default True.
    """
    assert(S.alphabet.letter_length == T.alphabet.letter_length)
    assert(S.alphabet.letters == T.alphabet.letters)
    letlen = S.alphabet.letter_length
    letters = S.alphabet.letters
    infostr, transcript = transcript.split(':', 1)
    transcript = transcript[1:] # skip the B
    indices, score = infostr.split('),', 1)
    idx_S, idx_T = indices[1:].split(',', 1)
    idx_S, idx_T = int(idx_S), int(idx_T)

    slines = tlines = []
    sline = tline = ''

    def print_lines(sline, tline, f):
        maxlen = max(len(sline), len(tline))
        sline, tline = sline.rjust(maxlen), tline.rjust(maxlen)
        f.write('%s\n%s\n' % (sline,tline))

    def new_line(sline, tline, _idx_S, _idx_T, f):
        print_lines(sline, tline, f)
        sline, tline = 'S[%d]: ' % _idx_S, 'T[%d]: ' % _idx_S
        return (max(len(sline), len(tline)), sline, tline)

    # The pre margin:
    pre_margin = min(margin, max(idx_S, idx_T))
    sline = 'S[%d]: ' % max(0, idx_S - pre_margin + 1)
    tline = 'T[%d]: ' % max(0, idx_T - pre_margin + 1)
    counter = max(len(sline), len(tline))
    for i in reversed(range(1, pre_margin)):
        if counter >= width:
            counter, sline, tline = new_line(sline, tline, idx_S+i, idx_T+i, f)
        sline += letters[S.c_idxseq[idx_S-i]] if i <= idx_S else ' '
        tline += letters[T.c_idxseq[idx_T-i]] if i <= idx_T else ' '
        counter += letlen

    # The alignment itself:
    for i,op in enumerate(transcript):
        if counter >= width:
            counter, sline, tline = new_line(sline, tline, idx_S, idx_T, f)
        if op in 'MS':
            s, t = letters[S.c_idxseq[idx_S]], letters[T.c_idxseq[idx_T]]
            idx_S += 1
            idx_T += 1
        elif op == 'I':
            s, t = '-', letters[T.c_idxseq[idx_T]]
            idx_T += 1
        elif op == 'D':
            s, t = letters[S.c_idxseq[idx_S]], '-'
            idx_S += 1
        on_color = color = None
        if colors:
            if op in 'MS':
                color = 'green' if op == 'M' else 'red'
            elif op == 'I':
                on_color = 'on_red'
            elif op == 'D':
                on_color = 'on_red'
        sline += colored(s, color=color, on_color=on_color)
        tline += colored(t, color=color, on_color=on_color)
        counter += letlen

    # The post margin:
    post_margin = min(margin, max(S.length-idx_S, T.length-idx_T))
    for i in range(post_margin):
        if counter >= width:
            counter, sline, tline = new_line(
                sline, tline, idx_S+i, idx_T+i, f)
        sline += letters[S.c_idxseq[idx_S+i]] if idx_S + i < S.length else ' '
        tline += letters[T.c_idxseq[idx_T+i]] if idx_T + i < T.length else ' '
        counter += letlen

    print_lines(sline, tline, f)
