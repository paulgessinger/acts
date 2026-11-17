import sys, collections


def load(p):
    c = collections.Counter()
    for l in open(p):
        c[tuple(l.split())] += 1
    return c


a, b = load(sys.argv[1]), load(sys.argv[2])
gained = sum((b - a).values())
lost = sum((a - b).values())
print(f"hits {sum(a.values())} -> {sum(b.values())}: gained {gained}, lost {lost}")
