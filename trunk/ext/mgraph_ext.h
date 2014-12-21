/*
 * File:   mgraph_ext.h
 * Author: bharath
 *
 * Created on 1 March, 2013, 10:18 PM
 */

#ifndef MGRAPH_EXT_H
#define	MGRAPH_EXT_H
#include "mgraph.h"
/*
    To use any igraph function
    mgraph mg;
    mg.(Call function)
    mg.createTopology_PastryHypercubeRing(int numnodes, int b, int numbits)
 */

void mgraph2graph(mgraph& mg, CDLib::graph& g) {
    g.clear();
    for (int k = 0; k < mg.numnodes(); k++)g.add_node();
    for (int k = 0; k < mg.numnodes(); k++)
        for (map<int, linkattribs>::iterator mi = mg.linkadjlist[k].begin(); mi != mg.linkadjlist[k].end(); mi++)
            g.add_edge(k, mi->first, 1);
}

void generate_mgraph_kademlia(graph& g, int nodes, int b) {
    mgraph mg;
    int bits = 128; //log2(nodes)+2;
    mg.createTopology_PastryHypercubeRing(nodes, b, bits);
    mgraph2graph(mg, g);
    g.set_graph_name("km_" + T2str<int>(nodes) + "_" + T2str<int>(b));
}

void generate_mgraph_kademlia_with_visibility(graph& g, int nodes, int b, double visibility) {
    generate_mgraph_kademlia(g, nodes, b);
    double removal = 1 - visibility;
    remove_edges_randomly(g, removal);
    string name = g.get_graph_name();
    g.set_graph_name(name + "_" + T2str<double>(removal));
}

#endif	/* MGRAPH_EXT_H */

