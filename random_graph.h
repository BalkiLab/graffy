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
#include "paths_and_components.h"
#include "random.h"

namespace CDLib {
    void generate_erdos_renyi_graph(graph& g, id_type num_nodes, double p);
    void generate_erdos_renyi_graph(graph& g, id_type num_nodes, id_type num_edges);
    //    Need re-implimentations
    void generate_scale_free_graph(graph& g, id_type num_nodes, id_type num_edges, double alpha, double beta);
    void generate_planted_partition_graph(graph& g, id_type num_comms, id_type comm_size, double pin, double pout, vector< node_set>& communities);
    long generate_lfr_graph(graph& g, id_type num_nodes, id_type num_edges, id_type max_degree, double tau, double tau2, double mixing_parameter, vector<node_set>& comms);
    //    Properly Implemented
    void generate_ring_graph(graph& g, id_type size);
    void generate_star_graph(graph& g, id_type size);
    void generate_clique_graph(graph& g, id_type size);
    void generate_spoke_graph(graph& g, id_type size);
    void generate_de_bruijn_graph(graph& g, id_type num_symbols, id_type sequence_length);
    void generate_chord_graph(graph& g, id_type num_nodes);
    void generate_LEET_chord_graph(graph& g, id_type num_nodes);
    void generate_kademlia_graph(graph& g, id_type num_nodes); // Doesn't work
    void generate_pref_attach_with_degree_seq(graph& g, vector<id_type>& deg_seq); // Need to be verified
    void init_empty_graph(graph& g, size_t size);
    //    Need to be verified
    void generate_configuration_model(graph& g, vector<id_type>& degree_sequence);
    void generate_prices_model(graph& g, size_t num_nodes, size_t num_of_out_degree, size_t in_degree_constant);
    void generate_barabasi_albert_model(graph& g, size_t num_nodes, size_t min_degree_of_node);
    bool generate_vertex_copying_model(graph& g, size_t num_nodes, size_t num_of_out_degree, size_t num_of_vertices_at_initial, double probability_to_copy_from_existing_vertex);
    bool generate_small_world_model(graph& g, size_t num_nodes, size_t degree_of_each_vertex, double probability_to_replace_edge);
    void generate_evolutionary_model_128_nodes_4_communities(vector<graph>& g, size_t num_inter_community_edges, size_t num_timesteps);
    void generate_evolutionary_model_800_nodes_4_communities(vector < vector < vector < double > > >& points, vector<double>& x_coordinates, vector<double>& y_coordinates, double variance, size_t num_timesteps);
    bool generate_ferrer_i_cancho_model(graph& g, size_t num_nodes, size_t max_failure_allowed, double degree_dist_controling_parameter, double probability_to_alter_edge, double initial_probability_of_edge);
    // Edges Rewiring
    id_type rewire_with_degree_distribution(graph& g, id_type degree_lower, id_type degree_higher, id_type max_trials, id_type num_rewires, int type);
    id_type improved_rewire_with_degree_distribution(graph& g, id_type degree_lower, id_type degree_higher, id_type max_trials, id_type num_rewires, int type);

    inline id_type rewire_with_degree_distribution(graph& g, id_type max_trials, id_type num_rewires) {
        return rewire_with_degree_distribution(g, 0, g.get_num_nodes(), max_trials, num_rewires, 0);
    }
    id_type random_rewire(graph& g, id_type degree_lower, id_type degree_higher, id_type max_trials, id_type num_rewires, int type);
    id_type random_rewire(graph& g, id_type num_rewires);
    id_type random_rewire(graph& g, double fraction_to_rewire);
    id_type randomize_chain_switching_method(graph& g, id_type max_trials, id_type num_rewires);

    inline id_type randomize_chain_switching_method(graph& g, double fraction) {
        if ((fraction < 0) && (fraction > 1) && (g.get_num_nodes() < 3) && (g.get_num_edges() < 2))
            return 0;
        return randomize_chain_switching_method(g, 100, (id_type) (fraction / g.get_num_edges()));
    }
    id_type make_quick_assortative(graph& g, id_type degree_cutoff, id_type max_trials, id_type num_rewires, int type);
    id_type make_quick_assortative(graph& g, id_type degree_cutoff, id_type num_rewires, bool random);
    id_type make_quick_disassortative(graph& g, id_type degree_cutoff, id_type num_rewires, bool random);
};

#endif	/* RANDOM_GRAPH_H */

