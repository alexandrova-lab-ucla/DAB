#!/usr/bin/env python3
# dijkstra's algorithm for finding paths from normalized covariance matrix of a protein simulation
# written nicely by wack Jack 20200222

import argparse


def main(csvfilename, res1, res2):
    dmat = []
    with open(csvfilename) as csvfile:
        for line in csvfile:
            line = line.split(',')
            newline = []
            for x in line:
                newline.append(1 - abs(float(x)))
            dmat.append(newline)
    paths = []
    distances = dmat[res1 - 1]
    set1 = set()
    set2 = set()
    for resnum in range(1, len(dmat) + 1):
        paths.append([res1, resnum])
        set2.add(resnum)
    while res2 not in set1:
        mindist = min([distances[x - 1] for x in set2])
        newres = distances.index(mindist) + 1
        for res in set2:
            newdist = dmat[newres - 1][res - 1] + mindist
            if newdist < distances[res - 1]:
                distances[res - 1] = newdist
                paths[res - 1] = paths[newres - 1] + [res]
        set2.remove(newres)
        set1.add(newres)
    path = paths[res2 - 1]
    pathcost = []
    for resind, res in enumerate(path[:-1]):
        pathcost.append(dmat[res - 1][path[resind + 1] - 1])
    return path, pathcost


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='analyzes shortest path between two residues using normalized '
                                                 'covariance matrix')

    parser.add_argument('csvfilename', type=str, help='name of csv file for normalized covariance matrix')
    parser.add_argument('res1', type=int, help='residue number (from DMD numbering scheme) for path start')
    parser.add_argument('res2', type=int, help='residue number (from DMD numbering scheme) for path end')

    args = parser.parse_args()

    print(main(args.csvfilename, args.res1, args.res2))
