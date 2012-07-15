/* 
 * File:   centrality.h
 * Author: bharath
 *
 * Created on 5 April, 2012, 2:36 AM
 */

#ifndef CENTRALITY_H
#define	CENTRALITY_H
#include "typedefs.h"
#include "graph.h"
#include "binary_heap.h"
using namespace std;
namespace CDLib 
{
    void betweeness_centralities(const graph& g, vector<double>& bc);
    double edge_clustering_coefficient(const graph&g,id_type from_id, id_type to_id);
    void degree_sequence(const graph& g, vector<id_type>& sequence);
    /*This gives the #nodes of degree indicated as the index of sequence variable in an undirected graph*/
};

#endif	/* CENTRALITY_H */

