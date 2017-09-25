import networkx as nx
import argparse
import matplotlib.pyplot as plt


def buildInitialNetwork(args):
    '''
    Build a networkx "Graph" object of the network
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
    f = plt.figure(figsize=(8, 8))
    a = f.add_subplot('111')
    colours = []
    for n in N.nodes(data=True):
        colours.append(colourdict[n[1]['host']])
    nx.draw_networkx(N, ax=a, node_color=colours)
    f.savefig("test.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--hostfile', dest='hostfile')
    parser.add_argument('--similarityfile', dest='similarityfile')
    parser.add_argument('--threshold', dest='threshold', type=float)
    parser.add_argument('--colourfile', dest='colours', type=str)
    args = parser.parse_args()
    colourdict = dict([line.strip().split("\t")
                       for line in open(args.colours).readlines()])
    N = buildInitialNetwork(args)
    plotInitialNetwork(N, args, colourdict)
