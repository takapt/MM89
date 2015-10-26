#! /usr/bin/env python3

import sys

def get_scores(filename):
    scores = {}
    for line in open(filename):
        s = line.split()
        seed = int(s[0])
        score = float(s[1])
        scores[seed] = score
    return scores

def get_seed_info(seed):
    try:
        f = open('inputs/' + str(seed))
    except:
        return ''

    f.readline()
    h = int(f.readline().strip())
    board = []
    for y in range(h):
        board.append(f.readline().strip())
    w = len(board[0])
    f.readline()
    max_changes = int(f.readline().strip())

    outside = sum(s.count('.') for s in board)
    inside = w * h - outside

    return '{:4d} {:3d} {:3d} {:3d} {:4d} {:6.3f} {:9d}'.format(w * h, w, h, max(w, h), max_changes, max_changes / inside, w * h * max_changes)

a, b = sys.argv[1], sys.argv[2]
p, q = get_scores(a), get_scores(b)
seeds = list(set(p) & set(q))
seeds = [seed for seed in seeds if seed > 0]
sum_ratio = 0

print('seed a b diff area w h max_wh changes changes_p areaxchange')
for seed in seeds:
    x, y = p[seed], q[seed]
    ratio = y - x

    sum_ratio += ratio

#     if ratio < 0.9:
#     if ratio > 1.2:
    print('{:>10} {:>10} {:>10} {:>7.3f}'.format(seed, x, y, ratio) + '      ' + get_seed_info(seed))

total_ratio = sum_ratio / len(seeds)
# print('total_ratio: {}'.format(total_ratio))

