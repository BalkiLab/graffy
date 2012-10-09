/* 
 * File:   random_graph.h
 * Author: bharath
 *
 * Created on 16 April, 2012, 4:25 PM
 */

#ifndef RANDOM_GRAPH_H
#define	RANDOM_GRAPH_H
#include "typedefs.h"
#include "graph.h"
//#include "extern.hpp"
#include "community_tools.h"

namespace CDLib {
    void generate_erdos_renyi_graph(graph& g, id_type num_nodes,double p);
    void generate_scale_free_graph(graph& g, id_type num_nodes,id_type num_edges,double alpha, double beta);
    void generate_planted_partition_graph(graph& g, id_type num_comms, id_type comm_size, double pin, double pout,vector< node_set>& communities);
    long generate_lfr_graph(graph& g,id_type num_nodes,id_type num_edges, id_type max_degree,double tau, double tau2,double mixing_parameter,vector<node_set>& comms);
    void generate_ring_graph(graph& g, id_type size);
    void generate_star_graph(graph& g, id_type size);
    void generate_clique_graph(graph& g, id_type size);
    void generate_spoke_graph(graph& g, id_type size);
    void generate_de_bruijn_graph(graph& g,id_type num_symbols,id_type sequence_length);
    void generate_chord_graph(graph& g,id_type num_nodes);
    void generate_LEET_chord_graph(graph& g,id_type num_nodes);
    void generate_kademlia_graph(graph& g,id_type num_nodes, id_type bucket_length);
    void generate_pref_attach_with_degree_seq(graph& g, vector<id_type>& deg_seq);
    void init_empty_graph(graph& g,size_t size);
    void generate_configuration_model(graph& g, vector<id_type>& degree_sequence);
    void generate_prices_model(graph& g,size_t num_nodes, size_t num_of_out_degree, size_t in_degree_constant);
    void generate_barabasi_albert_model(graph& g,size_t num_nodes,size_t min_degree_of_node);
    void generate_vertex_copying_model(graph& g,size_t num_nodes,size_t num_of_out_degree,size_t num_of_vertices_at_initial,double probability_to_copy_from_existing_vertex);
    void generate_small_world_model(graph& g,size_t num_nodes,size_t degree_of_each_vertex,double probability_to_replace_edge);
    void generate_evolutionary_model_128_nodes_4_communities(vector<graph>& g,size_t num_inter_community_edges,size_t num_timesteps);
    void generate_evolutionary_model_800_nodes_4_communities(vector < vector < vector < double > > >& points,vector<double>& x_coordinates,vector<double>& y_coordinates,double variance,size_t num_timesteps);
};

#endif	/* RANDOM_GRAPH_H */

