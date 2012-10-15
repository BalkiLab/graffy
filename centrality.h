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
#include "paths_and_components.h"

using namespace std;
namespace CDLib 
{
    void betweeness_centralities(const graph& g, vector<double>& bc);
    double edge_clustering_coefficient(const graph&g,id_type from_id, id_type to_id);
    void degree_centralities(const graph& g, vector<id_type>& degrees);
    void degree_centralities_normalized(const graph& g, vector<double>& degrees);
    void degree_sequence(const graph& g, vector<id_type>& sequence);
    /*This gives the #nodes of degree indicated as the index of sequence variable in an undirected graph*/
//    Overoaded Node Clustering Coefficient for single and all nodes.
    double node_clustering_coefficient(const graph&g, id_type node);
    void node_clustering_coefficient(const graph&g, vector<double> nodes);
    double closeness_centrality_original(const graph& g, id_type node);
    void closeness_centralities_original(const graph& g, vector<double>& closeness);
    double closeness_centrality(const graph& g, id_type node);
    void closeness_centralities(const graph& g, vector<double>& closeness);
    void eigenvector_centralities(const graph& g, vector<double>& eigenvector);
    void eigenvector_centralities_normalized(const graph& g, vector<double>& eigenvector);
};

#endif	/* CENTRALITY_H */

