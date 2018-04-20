#!/usr/bin/env
from math import sqrt
from matplotlib import pyplot as plt
import random
from itertools import product


def normalize(x, y):
    norm = sqrt(x**2 + y**2)
    return x / norm, y / norm


def rw(n_steps):
    x, y = 0, 0
    xs, ys = [], []
    for _ in range(n_steps):
        dx, dy = random.choice(moves_2d)
        x, y = x + dx, y + dy
        xs.append(x)
        ys.append(y)

    return xs, ys


if __name__ == '__main__':
    resolution = 20
    moves_1d = [i for i in range(-resolution, resolution) if i]
    moves_2d = [normalize(x, y) for x, y in product(moves_1d, moves_1d)]

    n_steps = 1000
    n_paths = 1000

    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(1, 1, 1)
    for i in range(n_paths):
        if (i + 1) % 100 == 0:
            print '%d/%d' % (i+1, n_paths)
        xs, ys = rw(n_steps)
        ax.plot(xs, ys, c='k', alpha=.1, lw=.5)

    ax.set_aspect('equal')
    ax.axis('off')
    fig.tight_layout()
    fig.savefig('blot.png', transparent=True, bbox_inches='tight',
                pad_inches=0.0, dpi=300)
