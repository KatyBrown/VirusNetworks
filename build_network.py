import networkx as nx
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os


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


def buildInitialNetwork(hostdict, scoredict):
    '''
    Build a networkx "Graph" object of the network.
    Nodes are coloured according to their host based on the
    file provided.
    '''
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


def plotNetwork(N, args, colourdict, outfile, lws=1):
    '''
    Draws a diagram of the network.
    '''
    f = plt.figure(figsize=(8, 8))
    a = f.add_subplot('111')
    colours = []
    for n in N.nodes(data=True):
        colours.append(colourdict[n[1]['host']])
    if args.shuffle:
        positions = shufflePositions(N, args.shuffle)
        nx.draw_networkx(N, ax=a, node_color=colours, pos=positions,
                         niter=args.niter, k=args.k, linewidths=lws)
    else:
        nx.draw_networkx(N, ax=a, node_color=colours,
                         niter=args.niter, k=args.k, linewidths=lws)
    a.collections[0].set_edgecolor("#000000") 
    a.axis('off')
    f.savefig(outfile)


def runConnectedComponents(N, args, colourdict, hostdict):
    '''
    Run analyses on the "connected component" subgraphs of the main graph.
    '''
    ccs = nx.connected_component_subgraphs(N)
    i = 1
    for cc in ccs:
        out = open("node_stats_%s.tsv" % i, "w")
        paths = nx.algorithms.shortest_paths.unweighted.all_pairs_shortest_path(N)
        dc = nx.degree_centrality(cc)
        bc = nx.betweenness_centrality(cc)

        m_dc = dc.keys()[0]
        m_bc = bc.keys()[0]
        for c in cc.nodes():
            if dc[c] > dc[m_dc]:
                m_dc = c
            if bc[c] > bc[m_bc]:
                m_bc = c
            out.write("%s\t%s\t%s\n" % (c, dc[c], bc[c]))
        out.close()
        lws = [1 if c != m_dc  else 3 for c in cc.nodes()]
        plotNetwork(cc, args, colourdict, "%s.png" % i, lws=lws)
        out = open("shortest_paths.tsv", "w")
        for path in paths:

            hosts = []
            for p in paths[path]:
                a_to_b = paths[path][p]
                hosts = set()
                for item in a_to_b:
                    h = hostdict[item]
                    hosts.add(h)
            
                hoststring = "->".join([hostdict[a] for a in a_to_b])

                out.write("%s\t%s\t%s\t%s\n" % (path, p, ",".join(a_to_b), hoststring))
        out.close()
        i += 1

    
def parseInput(args):
    '''
    Convert the input files into dictionaries.
    '''
    similarity = [line.strip().split("\t") for line in open(
        args.similarityfile).readlines()]
    scoredict = dict()
    for s in similarity:
        scoredict["%s_%s" % (s[0], s[1])] = float(s[2])

    hosts = [line.strip().split("\t") for line in open(
        args.hostfile).readlines()]
    hostdict = dict(hosts)

    colourdict = dict([line.strip().split("\t")
                       for line in open(args.colours).readlines()])

    return (hostdict, scoredict, colourdict)


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
    hostdict, scoredict, colourdict = parseInput(args)

    N = buildInitialNetwork(hostdict, scoredict)
    plotNetwork(N, args, colourdict, "network.png")
    runConnectedComponents(N, args, colourdict, hostdict)
