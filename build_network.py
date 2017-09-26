import networkx as nx
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import random
import matplotlib.patches as mpatches

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
    print ("Building initial network\n")
    N = nx.Graph()
    for s in scoredict:
        v1, v2 = s.split("|")
        if v1 not in N:
            N.add_node(v1, host=hostdict[v1])
        if v2 not in N:
            N.add_node(v2, host=hostdict[v2])
        N.add_edge(v1, v2, weight=scoredict[s])
    print ("Network built with %s nodes and %s edges" % (len(N.nodes()), len(N.edges())))
    return (N)


def plotNetwork(N, args, colourdict, outfile, lws=1):
    '''
    Draws a diagram of the network.
    '''
    print ("Plotting\n")
    f = plt.figure(figsize=(8, 8))
    a = f.add_subplot('111')
    colours = []
    hosts = []
    for n in N.nodes(data=True):
        colours.append(colourdict[n[1]['host']])
        hosts.append(n[1]['host'])
    hosts = set(hosts)
    patches = []
    labels = []
    for host in hosts:
        labels.append(host)
        patches.append(mpatches.Patch(color=colourdict[host]))
    if args.shuffle:
        positions = shufflePositions(N, args.shuffle)
        nx.draw_networkx(N, ax=a, node_color=colours, pos=positions,
                         niter=args.niter, k=args.k, linewidths=lws,
                         with_labels=False)
    
    else:
        nx.draw_networkx(N, ax=a, node_color=colours,
                         niter=args.niter, k=args.k, linewidths=lws,
                         with_labels=False)
    a.collections[0].set_edgecolor("#000000") 
    a.axis('off')
    a.legend(handles=patches, labels=labels)
    f.savefig(outfile)

def getDCBC(N):
    '''
    Calculated betweenness centrality and degree centrality
    '''
    dc = nx.degree_centrality(N)
    bc = nx.betweenness_centrality(N)
    m_dc = dc.keys()[0]
    m_bc = bc.keys()[0]
    for c in N.nodes():
        if dc[c] > dc[m_dc]:
            m_dc = c
        if bc[c] > bc[m_bc]:
            m_bc = c
    
    m_dc = dc.keys()[0]
    m_bc = bc.keys()[0]
    return (dc, bc, m_dc, m_bc)

def writeNodeStats(N, dc, bc, outfile):
    out = open(outfile, "w")
    out.write("Node\tDegree_Centrality\tBetweenness_Centrality\n")
    for c in N.nodes():
        out.write("%s\t%s\t%s\n" % (c, dc[c], bc[c]))    
    out.close()

def runConnectedComponents(N, args, colourdict, hostdict):
    '''
    Run analyses on the "connected component" subgraphs of the main graph.
    '''
    ccs = nx.connected_component_subgraphs(N)
    i = 1
    print ("Analysing %s connected components\n" % len(list(ccs)))
    ccs = nx.connected_component_subgraphs(N)
    for cc in ccs:
        print ("Plotting connected component %i" % i)

        
        dc, bc, m_dc, m_bc = getDCBC(cc)
        writeNodeStats(cc, dc, bc, "node_%i.tsv" % i)
        lws = [1 if c != m_dc  else 3 for c in cc.nodes()]
        plotNetwork(cc, args, colourdict, "%s.png" % i, lws=lws)

        for c in cc.nodes():
            if c == m_dc or c == m_bc:
                seq = fastadict[c]
                out = open("temp.fasta", "w")
                out.write(">c\n%s\n" % seq)
                out.close()
                os.system('../usearch10.0.240_i86linux32 -ublast temp.fasta -db %s -evalue 1e-9 -strand both -blast6out temp_%i.out' % (args.ublast, i))
        plotNetwork(cc, args, colourdict, "%s.png" % i, lws=lws)
        paths = nx.algorithms.shortest_paths.unweighted.all_pairs_shortest_path(N)
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

def getColours():
    return (['#e6194b', '#3cb44b', '#ffe119', '#0082c8',
             '#f58231', '#911eb4', '#46f0f0', '#f032e6',
             '#d2f53c', '#fabebe', '#008080', '#e6beff',
             '#aa6e28', '#fffac8', '#800000', '#aaffc3',
             '#808000', '#ffd8b1', '#000080', '#808080',
             '#FFFFFF', '#000000'])

def randomColour():
    colour = '#{:02x}{:02x}{:02x}'.format(*map(
        lambda x: random.randint(0, 255), range(3)))
    return (colour)

def parseInput(args):
    '''
    Convert the input files into dictionaries.
    '''
    print ("Importing distance matrix\n")
    scoredict = dict()
    hostdict = dict()
    familydict = dict([line.split("\t")[0:2] for line in open(args.hostfile).readlines()])
    with open(args.similarityfile) as similarity:
        for s in similarity:
            s = s.split("\t")
            if float(s[2]) >= args.threshold:
                host1 = s[0][0:4]
                host2 = s[1][0:4]
                family1 = familydict[host1]
                family2 = familydict[host2]
                one = (family1 != family2 and args.transtype == "crossfamily")
                two = (host1 != host2 and args.transtype == "crosshost")
                three = (args.transtype == "all")
                if (one or two or three):
                    scoredict["%s|%s" % (s[0], s[1])] = float(s[2])
                    hostdict[s[0]] = family1
                    hostdict[s[1]] = family2
    print ("Imported %i datapoints from distance matrix" % len(scoredict))
    colours = getColours()
    colourdict = dict()
    i = 0
    for host in familydict:
        colourdict[familydict[host]] = randomColour()
        i += 1
    return (hostdict, scoredict, colourdict)

def getFastaDict(fasta):
    print ("Parsing FASTA file to dictionary")
    fastadict = dict()
    seq = []
    nam = 0
    with open(fasta) as infile:
        for line in infile:
            line = line.strip()
            if line[0] == ">":
                if nam != 0:
                    fastadict[nam] = "".join(seq)
                    seq = []
                nam = line.replace(">", "")
            else:
                seq.append(line)
    fastadict[nam] = "".join(seq)
    print ("Generated FASTA dictionary of length %i" % len(fastadict))
    return (fastadict)


if __name__ == "__main__":
    random.seed(10)
    parser = argparse.ArgumentParser()

    parser.add_argument('--hostfile', dest='hostfile', type=str)
    parser.add_argument('--similarityfile', dest='similarityfile',
                        type=str)
    parser.add_argument('--fasta', dest='fasta', type=str)
    parser.add_argument('--threshold', dest='threshold', type=float)
    parser.add_argument('--ublast', dest='ublast', type=str)

    parser.add_argument('--shuffle', dest='shuffle', type=float)
    parser.add_argument('--nodedist', dest='k', type=float,
                        default=0.1)
    parser.add_argument('--niter', dest='niter', type=int,
                        default=1000)
    parser.add_argument('--transtype', dest='transtype', type=str,
                        default='crossfamily')

    args = parser.parse_args()
    hostdict, scoredict, colourdict = parseInput(args)
    fastadict = getFastaDict(args.fasta)
    N = buildInitialNetwork(hostdict, scoredict)
    plotNetwork(N, args, colourdict, "network.png")
    runConnectedComponents(N, args, colourdict, hostdict)
