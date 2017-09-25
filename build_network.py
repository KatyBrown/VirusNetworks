import networkx as nx
import argparse


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
    print (N.nodes())

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--hostfile', dest='hostfile')
    parser.add_argument('--similarityfile', dest='similarityfile')
    parser.add_argument('--threshold', dest='threshold', type=float)
    args = parser.parse_args()

    buildInitialNetwork(args)
