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


def plotNetwork(N, args, colourdict, outfile, lws=1, txt=None, title=None):
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
    if txt is not None:
        a.text(0, -0.2, txt)
    if title is not None:
        a.set_title(title)
    a.legend(handles=patches, labels=labels, loc=1)
    f.savefig(outfile)
    f.savefig(outfile.replace(".png", ".svg"))

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
    return (dc, bc, m_dc, m_bc)


def runUBLAST(args, ID, fastadict):
    seq = fastadict[ID]
    out = open("temp.fasta", "w")
    out.write(">c\n%s\n" % seq)
    out.close()
    os.system('%s -ublast temp.fasta -db %s -evalue 1E-9 \
    -strand both -blast6out temp.out -qmask none -quiet' % (args.path_ublast, args.ublast))
    tab = open("temp.out").readline().split("\t")
    if len(tab) != 1:
        string = "%s_%s_%.2f" % (tab[1], tab[3], float(tab[2]))
    else:
        string = 0
    os.unlink("temp.out")
    os.unlink("temp.fasta")
    return (string)

def outputFasta(N, fastadict, outfile):
    out = open(outfile, "w")
    for n in N.nodes():
        out.write(">%s\n%s\n" % (n, fastadict[n]))
    out.close()

def runConnectedComponents(N, args, colourdict, hostdict):
    '''
    Run analyses on the "connected component" subgraphs of the main graph.
    '''
    ccs = nx.connected_component_subgraphs(N)
    i = 1
    print ("Analysing %s connected components\n" % len(list(ccs)))
    ccs = nx.connected_component_subgraphs(N)
    blast_results = []
    shortest_paths = []
    centrality = []
    for cc in ccs:
        
        dc, bc, m_dc, m_bc = getDCBC(cc)
        lws = []
        for c in cc.nodes():
            nodestats = "%i\t%s\t%.3f\t%.3f\n" % (i, c, dc[c], bc[c])
            centrality.append(nodestats)

            if c == m_bc:
                dc_string = runUBLAST(args, c, fastadict)
            if c == m_bc:
                bc_string = runUBLAST(args, c, fastadict)
            if c == m_dc or c == m_bc:
                lws.append(3)
            else:
                lws.append(1)
        if dc_string == 0 and bc_string == 0:
            for c in cc.nodes():
                other_string = runUBLAST(args, c, fastadict)
                if other_string != 0:
                    break
        else:
            other_string = 0
        if (bc_string != 0 or dc_string != 0
            or other_string != 0) and len(cc.nodes()) > 2:
            bc_string, dc_string, other_string = str(bc_string), str(
                dc_string), str(other_string)
            blast_results.append("%i\t%s\t%s\t%s\n" % (i, bc_string, dc_string, other_string))
            blast_result = "\n".join([bc_string, dc_string, other_string])
            print ("Plotting connected component %i" % i)
            plotNetwork(cc, args, colourdict,
                        "%s.png" % i, lws=lws,
                        txt=blast_result, title="node_%i" % i)
            outputFasta(cc, fastadict, "%s.fasta" % i)
            paths = nx.algorithms.shortest_paths.unweighted.all_pairs_shortest_path(N)

            for path in paths:
                hosts = []
                for p in paths[path]:
                    a_to_b = paths[path][p]
                    hosts = set()
                    for item in a_to_b:
                        h = hostdict[item]
                        hosts.add(h)
            
                    hoststring = "->".join([hostdict[a] for a in a_to_b])
                    shortest_paths.append("%i\t%s\t%s\t%s\t%s\n" % (
                        i, path, p, ",".join(a_to_b), hoststring))
            i += 1
    out = open("ublast_results.tsv", "w")
    out.write("connected_component\tbest_hit_betweenness\tbest_hit_degree\
    \tbest_hit_other\n")
    for line in blast_results:
        out.write(line)
    out.close()
    out = open("node_centrality.tsv", "w")
    for line in centrality:
        out.write(line)
    out.close()
    out = open("shortest_paths.tsv", "w")
    for line in shortest_paths:
        out.write(line)
    out.close()

def getColours():
    return (["#6A3A4C", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46",
             "#008941", "#006FA6", "#A30059", "#FFDBE5", "#7A4900",
             "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF",
             "#997D87", "#5A0007", "#809693", "#FEFFE6", "#1B4400",
             "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80", "#61615A",
             "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9",
             "#B903AA", "#D16100", "#DDEFFF", "#000035", "#7B4F4B",
             "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
             "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744",
             "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062",
             "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1",
             "#788D66", "#885578", "#FAD09F", "#FF8A9A", "#D157A0",
             "#BEC459", "#456648", "#0086ED", "#886F4C", "#34362D",
             "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9",
             "#FF913F", "#938A81", "#575329", "#00FECF", "#B05B6F",
             "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
             "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C",
             "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F",
             "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465",
             "#922329", "#5B4534", "#FDE8DC", "#404E55", "#0089A3",
             "#CB7E98", "#A4E804", "#324E72"])


def randomColour():
    colours = getColours()
    ind = np.random.randint(0, 64)
    return (colours[ind])

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
    np.random.seed(16)
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
    parser.add_argument("--ublastpath", dest='path_ublast',
                        type=str, default='ublast')

    args = parser.parse_args()
    hostdict, scoredict, colourdict = parseInput(args)
    fastadict = getFastaDict(args.fasta)
    N = buildInitialNetwork(hostdict, scoredict)
    plotNetwork(N, args, colourdict, "network.png", title="Full Network")
    runConnectedComponents(N, args, colourdict, hostdict)
