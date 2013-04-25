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
#include "datastructures.h"
#include "paths_and_components.h"

using namespace std;
namespace CDLib {
    void betweeness_centralities(const graph& g, vector<double>& bc);
    void betweenness_centralities_normalized(const graph& g, vector<double>& bcn);
    double edge_clustering_coefficient(const graph&g, id_type from_id, id_type to_id);

    template <typename T>
    void degree_centralities(const graph& g, vector<T>& degrees) {
        //    Return the out-degree centrality of all the nodes in the given graph.
        degrees.clear();
        for (id_type i = 0; i < g.get_num_nodes(); i++)
            degrees.push_back((T) g.get_node_out_degree(i));
    }
    void degree_centralities_normalized(const graph& g, vector<double>& degrees);
    void degree_sequence(const graph& g, vector<id_type>& sequence);
    /*This gives the #nodes of degree indicated as the index of sequence variable in an undirected graph*/
    //    Overoaded Node Clustering Coefficient for single and all nodes.
    double node_clustering_coefficient(const graph&g, id_type node);
    void node_clustering_coefficient(const graph&g, vector<double>& nodes);
    void node_clustering_coefficient_normalized(const graph& g, vector<double>& nodes);
    double average_clustering_coefficient(const graph& g);
    double closeness_centrality_original(const graph& g, id_type node);
    void closeness_centralities_original(const graph& g, vector<double>& closeness);
    double closeness_centrality(const graph& g, id_type node);
    void closeness_centralities(const graph& g, vector<double>& closeness);
    void closeness_centralities_normalized(const graph& g, vector<double>& closeness);
    void eigenvector_centralities(const graph& g, vector<double>& eigenvector);
    void eigenvector_centralities_normalized(const graph& g, vector<double>& eigenvector);
    double efficiency_centrality(const graph& g, id_type node);
    void efficiency_centralities(const graph& g, vector<double>& centralities);
    // Below functions return <node label, value> pair
    pair<string, double> get_max_degree_node(const graph & g);
    pair<string, double> get_max_degree_node(const graph & g, const node_set_string& elements);
    pair<string, double> get_max_betweenness_node(const graph & g);
    pair<string, double> get_max_betweenness_node(const graph & g, const node_set_string& elements);
    pair<string, double> get_max_efficiency_centrality_node(const graph & g);
    pair<string, double> get_max_efficiency_centrality_node(const graph & g, const node_set_string& elements);
};

#endif	/* CENTRALITY_H */

