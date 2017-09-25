import networkx as nx
import argparse
import matplotlib.pyplot as plt
import numpy as np


def shufflePositions(N, level):
    '''
    Shuffles the starting positions in the diagram of the network - used to
    tweak the layout so all the nodes are visible.
    Level is a float specifying how far apart to put the initial points.
    The networkx spring layout function is then applied.
    '''
    cliques = list(nx.find_cliques(N))
    positions = dict()
    r = range(len(cliques))
    xpos = np.random.choice(r, len(cliques), replace=False)
    ypos = np.random.choice(r, len(cliques), replace=False)
    i = 0
    for L in list(nx.find_cliques(N)):
        xp = xpos[i]
        yp = ypos[i]
        x_range = np.linspace(xp, xp + level, len(L))
        y_range = np.linspace(yp, yp + level, len(L))
        j = 0
        for l in L:
            positions[l] = [x_range[j], y_range[j]]
            j += 1
        i += 1
    return positions


def buildInitialNetwork(args):
    '''
    Build a networkx "Graph" object of the network.
    Nodes are coloured according to their host based on the
    file provided.
    '''
    hosts = [line.strip().split("\t") for line in open(
        args.hostfile).readlines()]
    hostdict = dict(hosts)
    similarity = [line.strip().split("\t") for line in open(
        args.similarityfile).readlines()]
    scoredict = dict()
    for s in similarity:
        scoredict["%s_%s" % (s[0], s[1])] = float(s[2])
    N = nx.Graph()
    for s in scoredict:
        if scoredict[s] >= args.threshold:
            v1, v2 = s.split("_")
            if v1 not in N:
                N.add_node(v1, host=hostdict[v1])
            if v2 not in N:
                N.add_node(v2, host=hostdict[v2])
            N.add_edge(v1, v2, weight=scoredict[s])
    return (N)


def plotInitialNetwork(N, args, colourdict):
    '''
    Draws an initial diagram of the whole network.
    '''
    f = plt.figure(figsize=(8, 8))
    a = f.add_subplot('111')
    colours = []
    for n in N.nodes(data=True):
        colours.append(colourdict[n[1]['host']])
    if args.shuffle:
        positions = shufflePositions(N, args.shuffle)
        nx.draw_networkx(N, ax=a, node_color=colours, pos=positions,
                         niter=args.niter, k=args.k)
    else:
        nx.draw_networkx(N, ax=a, node_color=colours,
                         niter=args.niter, k=args.k)
    a.axis('off')
    f.savefig("test.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--hostfile', dest='hostfile', type=str)
    parser.add_argument('--similarityfile', dest='similarityfile',
                        type=str)
    parser.add_argument('--threshold', dest='threshold', type=float)
    parser.add_argument('--colourfile', dest='colours', type=str)
    parser.add_argument('--shuffle', dest='shuffle', type=float)
    parser.add_argument('--nodedist', dest='k', type=float,
                        default=0.1)
    parser.add_argument('--niter', dest='niter', type=int,
                        default=1000)
    args = parser.parse_args()
    colourdict = dict([line.strip().split("\t")
                       for line in open(args.colours).readlines()])
    N = buildInitialNetwork(args)
    plotInitialNetwork(N, args, colourdict)
