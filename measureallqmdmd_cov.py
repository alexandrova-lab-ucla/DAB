#!/usr/bin/env python3
# makes covariance matrix for everyone
# written by Jack 20190122 edited 20200521

import argparse
from contextlib import contextmanager
import csv
import numpy as np
import os
import os.path
import subprocess

# only path to change (probably)
dmdpath = '/u/project/ana/ana/bin/DMD'  # for old DMD
# dmdpath = '/u/home/j/jtfuller/jackscripts/pdmd'  # for Hoffmann
if not os.path.isdir(dmdpath):  # actually Ima make this automatically change the path for XSEDE
    dmdpath = '/home/jtfuller/jackscripts/pdmd'  # for XSEDE

dmdparpath = dmdpath + '/parameter'
dmdbinpath = dmdpath + '/bin'  # where complex.linux, pdmd.linux, and complex_M2P.linux live


@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)


def main(jobpath):
    jobpath = jobpath.rstrip('/')
    jobname = jobpath.split('/')[-1]
    repname = jobname
    repnum = 1
    snapshotnum = 0

    # loops through all replicates
    while os.path.isdir(jobpath + "/" + repname):
        repdir = '{0}/{1}/{1}/'.format(jobpath, repname)
        iteration = 0
        missingiter = 0
        # loops through all iterations of the replicate until three in a row are missing
        while missingiter < 3:
            iterdir = '{}Iteration_{}/'.format(repdir, iteration)
            if os.path.exists(iterdir):
                missingiter = 0
                # makes movie.pdb
                with cd(iterdir):
                    subprocess.run([os.path.join(dmdbinpath, 'complex_M2P.linux'), dmdparpath, 'new.pdb',
                                    'topparam.list', 'movie', 'movie.pdb', 'inConstr'], stdout=subprocess.DEVNULL,
                                   stderr=subprocess.DEVNULL)
                moviepath = '{}movie.pdb'.format(iterdir)
                if os.path.exists(moviepath):
                    pdblines = []
                    with open(moviepath) as moviefile:
                        # loops through all snapshots of the iteration
                        for line in moviefile:
                            if not line.startswith('ENDMDL'):
                                # collects lines until 'ENDMDL'
                                pdblines.append(line)
                            else:
                                snapshotnum += 1
                                alphacoords = []
                                for pdbline in pdblines:
                                    if pdbline[12:16].strip() == 'CA':
                                        coords = np.asarray(
                                            [float(pdbline[30:38]), float(pdbline[38:46]), float(pdbline[46:54])])
                                        coords -= np.asarray([100, 100, 100])
                                        alphacoords.append(coords)
                                if snapshotnum == 1:
                                    rirj = np.zeros([len(alphacoords), len(alphacoords)])
                                    ri = np.zeros([len(alphacoords), 3])
                                    ri2 = np.zeros(len(alphacoords))
                                for idx, alphacoord in enumerate(alphacoords):
                                    for idx2, alphacoord2 in enumerate(alphacoords):
                                        if idx2 >= idx:
                                            rirj[idx][idx2] += np.dot(alphacoord, alphacoord2)
                                    ri[idx] += alphacoord
                                    ri2[idx] += np.dot(alphacoord, alphacoord)
                                # reset pdblines for next snapshot
                                pdblines = []
                    os.remove(moviepath)
                else:
                    print(' movie.pdb not found for {}_{}'.format(repname, iteration))
            else:
                missingiter += 1
            iteration += 1
        if repname == jobname:
            repname = jobname + "2"
            repnum = 2
        else:
            repnum += 1
            repname = jobname + str(repnum)

    # Now average everything
    rirj = rirj / snapshotnum
    ri = ri / snapshotnum
    ri2 = ri2 / snapshotnum

    # Now compute covariance matrix
    cov = np.zeros([len(rirj), len(rirj)])
    for idx in range(len(rirj)):
        for idx2 in range(len(rirj)):
            if idx2 >= idx:
                cov[idx][idx2] = (rirj[idx][idx2] - np.dot(ri[idx], ri[idx2]))\
                                 / np.sqrt((ri2[idx] - np.dot(ri[idx], ri[idx]))
                                           * (ri2[idx2] - np.dot(ri[idx2], ri[idx2])))
            else:
                cov[idx][idx2] = cov[idx2][idx]

    # Now print to csv
    with open(jobname + '_cov.csv', 'w') as csvfile:
        writer = csv.writer(csvfile)
        for row in cov:
            writer.writerow(row)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='calculates covariance matrix for science')

    parser.add_argument('jobpath', help='directory where all replicates are found')

    args = parser.parse_args()

    main(args.jobpath)
