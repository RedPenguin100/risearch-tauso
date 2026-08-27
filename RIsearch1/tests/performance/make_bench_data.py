#!/usr/bin/env python3
"""Write the FASTA files the README's performance table is measured on.

Deterministic: the seed is fixed, so this produces the same bytes anywhere and
the published timings can be reproduced.

Every target file holds the same 2 400 000 nt and differs only in how it is
divided into records, so the alignment work is the same in each. Every query is
a 20 nt reverse complement of a stretch of that sequence, so the queries pair
against the target rather than scoring like noise.

    python3 make_bench_data.py
    RIsearch -q q_16.fa -t t_2500.fa -p3 > /dev/null
"""
import random

TOTAL_NT = 2_400_000
RECORD_LENGTHS = (600, 2500, 20000, TOTAL_NT)
QUERY_COUNTS = (16, 64, 256, 2000)
QUERY_NT = 20
COMPLEMENT = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}

random.seed(1)
sequence = ''.join(random.choice('ACGU') for _ in range(TOTAL_NT))

for record_length in RECORD_LENGTHS:
    with open('t_%d.fa' % record_length, 'w') as out:
        for i in range(TOTAL_NT // record_length):
            out.write('>r%06d\n%s\n'
                      % (i, sequence[i * record_length:(i + 1) * record_length]))

for count in QUERY_COUNTS:
    with open('q_%d.fa' % count, 'w') as out:
        for i in range(count):
            # A prime stride so the queries are spread over the whole sequence.
            start = (i * 7919) % (TOTAL_NT - QUERY_NT)
            site = sequence[start:start + QUERY_NT]
            out.write('>aso%05d\n%s\n'
                      % (i, ''.join(COMPLEMENT[c] for c in reversed(site))))
