/* 
 * File:   divisive_algorithms.h
 * Author: bharath
 *
 * Created on 5 April, 2012, 12:42 PM
 */

#ifndef DIVISIVE_ALGORITHMS_H
#define	DIVISIVE_ALGORITHMS_H
#include "graph.h"
#include "centrality.h"
#include "paths_and_components.h"
#include "typedefs.h"
using namespace std;

namespace CDLib 
{
    typedef vector< vector<node_set> > dendrogram;
    struct edge_radicchi
    {
        id_type from_id;
        id_type to_id;
        double ecc;
        edge_radicchi();
        edge_radicchi(id_type from,id_type to,double ec);
    };
    void girvan_newman_2002(const graph& g, dendrogram& dendro);
    void radicchi_et_al_2004(const graph& g, dendrogram & dendro);
    void vragovic_et_al_2006(const graph& g, dendrogram & dendro);
    void rattigan_et_al_2007(const graph& g, dendrogram & dendro);
};

#endif	/* DIVISIVE_ALGORITHMS_H */

