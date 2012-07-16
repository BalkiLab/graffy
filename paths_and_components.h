/* 
 * File:   community_tools.h
 * Author: bharath
 *
 * Created on April 2, 2012, 3:18 PM
 */

#ifndef PATHS_AND_COMPONENTS_H
#define	PATHS_AND_COMPONENTS_H

#include "graph.h"
#include "disjoint_set.h"
#include "binary_heap.h"
#include "random_graph.h"

using namespace std;

namespace CDLib
{
    void dfs_visitor(const graph& g, node_set& visited, id_type source);
    void bfs_visitor(const graph& g, node_set& visited, id_type source);
    
    bool is_path_present(const graph& g, id_type source, id_type dest);
    
    id_type get_component_around_node_undirected(const graph& g,id_type id, node_set& visited);
    id_type get_component_around_node_weak(const graph&g,id_type id,node_set& visited);
    id_type get_component_around_node_strong(const graph& g, id_type id, node_set& visited);

    bool is_connected_undirected(const graph& g);
    bool is_connected_weakly(const graph& g);
    bool is_connected_strongly(const graph& g);

    id_type get_connected_components_undirected(const graph& g, vector<node_set>& components);
    id_type get_weakly_connected_components(const graph& g, vector<node_set>& components);
    id_type get_strongly_connected_components(const graph& g, vector<node_set>& components);

    bool has_negative_edge_weights(const graph& g);

    double single_source_shortest_paths_bfs(const graph& g,id_type source,vector<double>& distances,vector< vector<id_type> >& preds);
    double single_source_shortest_paths_djikstra(const graph&g,id_type source,vector<double>& distances,vector< vector<id_type> >& preds);
    double diameter(const graph& g);
    void all_pairs_shortest_paths(const graph& g, vector< vector<double> >& path_matrix);
    void single_source_shortest_paths_djikstra_with_paths(const graph&g,id_type source,vector<double>& distances,vector< vector<id_type> >& paths);
    bool get_topological_ordering(const graph& g,vector<id_type>& ordering);
    
    void get_all_paths(const graph& g,id_type source, id_type dest,vector<id_type>& paths);
    void all_path_lenth_Monte_Carlo(const graph& g, vector< vector<double> >& paths, long monte_c);
    double blocking_probability(id_type number_of_nodes, id_type degree, id_type visited);
    void alternate_path_length_destabilization(graph&g,id_type source,vector<double>& alternate_distances);
    void generate_ferrer_i_cancho_model(graph& g,size_t num_nodes, size_t max_failure_allowed,double degree_dist_controling_parameter,double probability_to_alter_edge,double initial_probability_of_edge);
    double efficiency_sw_global(graph& g, bool type);
    double efficiency_sw_global_monte_carlo(graph& g);
    
    double connectivity_entropy(graph& g);
    double path_entropy(graph& g);
    
    id_type hop_distance_matrix(const graph& g, vector< vector<id_type> > & path_mat);
    


};

#endif	/* COMMUNITY_TOOLS_H */

