/*
 * File:   igraph_ext.h
 * Author: bharath
 *
 * Created on 1 March, 2013, 10:16 PM
 */

#ifndef IGRAPH_EXT_H
#define	IGRAPH_EXT_H

#include "includes.h"
#include "igraph.h"

/*
        To use any igraph function
        igraph_ext obj;
        obj.(Call functions);
 */


class igraph_ext {
private:
    igraph_t ig;
    bool ig_status;
    igraph_matrix_t merges;
    igraph_vector_t membership;
    igraph_vector_t modularity;
    igraph_arpack_options_t options;
    igraph_real_t codelength;
    string graph_name;
public:

    igraph_ext() {
        ig_status = false;
        init_igraph_arrays();
    }

    int import_graph(const graph& g) {
        igraph_vector_t edges = IGRAPH_VECTOR_NULL;
        IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
        IGRAPH_CHECK(igraph_vector_reserve(&edges, g.get_num_edges()));
        for (id_type i = 0; i < g.get_num_nodes(); i++) {
            for (adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i); aeit++) {
                if (i <= aeit->first) {
                    IGRAPH_CHECK(igraph_vector_push_back(&edges, static_cast<long> (i)));
                    IGRAPH_CHECK(igraph_vector_push_back(&edges, static_cast<long> (aeit->first)));
                }
            }
        }
        IGRAPH_CHECK(igraph_create(&ig, &edges, g.get_num_nodes(), g.is_directed()));
        graph_name = g.get_graph_name();
        ig_status = true;
        igraph_vector_destroy(&edges);
        IGRAPH_FINALLY_CLEAN(1);
        return 1;
    }

    void init_igraph_arrays() {
        igraph_vector_init(&membership, 0);
        igraph_vector_init(&modularity, 0);
        igraph_matrix_init(&merges, 0, 0);
        igraph_arpack_options_init(&options);
    }

    igraph_ext(const graph& g) {
        import_graph(g);
        init_igraph_arrays();
    }

    void return_community_structure(vector<node_set>& comms) {
        comms.clear();
        vector<id_type> labels(igraph_vcount(&ig), 0);
        for (long i = 0; i < igraph_vector_size(&membership); i++)labels[i] = VECTOR(*&membership)[i];
        convert_labels_to_communities(labels, comms);
    }

    pair<long, long> graph_details() {
        return make_pair(igraph_vcount(&ig), igraph_ecount(&ig));
    }

    int rewire_edges(double probability) {
        int ig_code = igraph_rewire_edges(&ig, probability, 0, 0);
        if (ig_code == 0) {
            graph_name = graph_name + "_ire" + T2str<double>(probability);
        }
        return ig_code;
    }

    int rewire_with_degree_distribution(int trials) {
        int ig_code = igraph_rewire(&ig, trials, IGRAPH_REWIRING_SIMPLE);
        //        if (ig_code == 0) {
        //            graph_name = graph_name + "_irdd" + T2str<int>(trials);
        //        }
        return ig_code;
    }

    int generate_configuration_model(vector<id_type>& degrees_sequence) {
        if (ig_status) igraph_destroy(&ig);
        igraph_vector_t igraph_degrees;
        igraph_vector_init(&igraph_degrees, degrees_sequence.size());
        for (id_type i = 0; i < degrees_sequence.size(); i++) VECTOR(igraph_degrees)[i] = static_cast<igraph_real_t> (degrees_sequence[i]);
        int error_code = igraph_degree_sequence_game(&ig, &igraph_degrees, 0, IGRAPH_DEGSEQ_VL);
        if (!error_code) {
            graph_name = "degree_sequence_" + T2str<id_type > (degrees_sequence.size());
            ig_status = true;
        }
        igraph_vector_destroy(&igraph_degrees);
        return error_code;
    }

    int generate_scale_free(id_type nodes, id_type edges, double exp) {
        if (ig_status) igraph_destroy(&ig);
        int error_code = igraph_barabasi_game(&ig, nodes, exp, edges, 0, 1, 1, 0, IGRAPH_BARABASI_PSUMTREE, 0);
        if (!error_code) {
            graph_name = "scale_free_" + T2str<id_type > (nodes) + "_" + T2str<id_type > (edges) + "_" + T2str<double>(exp);
            ig_status = true;
        }
        return error_code;
    }

    int generate_small_world(id_type size, id_type dimension, double probability) {
        if (ig_status) igraph_destroy(&ig);
        int error_code = igraph_watts_strogatz_game(&ig, dimension, size, 1, probability, false, false);
        if (!error_code) {
            graph_name = "sw_" + T2str<id_type > (size) + "_" + T2str<id_type > (dimension) + "_" + T2str<double>(probability);
            ig_status = true;
        }
        return error_code;
    }

    int generate_forest_fire_model(id_type nodes, double forward_probability, double backward_ratio, id_type pambs, bool directed) {
        if (ig_status) igraph_destroy(&ig);
        int error_code = igraph_forest_fire_game(&ig, nodes, forward_probability, backward_ratio, pambs, directed);
        if (!error_code) {
            graph_name = "ff_" + T2str<id_type > (nodes) + "_" + T2str<double> (forward_probability) + "_" + T2str<double>(backward_ratio);
            ig_status = true;
        }
        return error_code;
    }

    int generate_de_brujin(id_type nodes, id_type symbols, id_type seq_ln) {
        if (ig_status) igraph_destroy(&ig);
        int error_code = igraph_de_bruijn(&ig, symbols, seq_ln);
        if (!error_code) {
            graph_name = "db_" + T2str<id_type > (symbols) + "_" + T2str<id_type > (seq_ln);
            ig_status = true;
        }
        return error_code;
    }

    int generate_er_graph(id_type nodes, id_type edges) {
        if (ig_status) igraph_destroy(&ig);
        int error_code = igraph_erdos_renyi_game(&ig, IGRAPH_ERDOS_RENYI_GNM, nodes, (double) edges, 0, 0);
        if (!error_code) {
            graph_name = "er_" + T2str<id_type > (nodes) + "_" + T2str<id_type > (edges);
            ig_status = true;
        }
        return error_code;
    }

    int generate_er_graph(id_type nodes, double probability) {
        if (ig_status) igraph_destroy(&ig);
        int error_code = igraph_erdos_renyi_game(&ig, IGRAPH_ERDOS_RENYI_GNP, nodes, probability, 0, 0);
        if (!error_code) {
            graph_name = "er_" + T2str<id_type > (nodes) + "_" + T2str<double>(probability);
            ig_status = true;
        }
        return error_code;
    }

    void leading_eigenvector_newman() {
        igraph_community_leading_eigenvector(&ig, &merges, &membership, igraph_vcount(&ig), &options, NULL, 0, NULL, NULL, NULL, NULL, NULL);
    }

    void infomap_cda(int num_trials) {
        igraph_community_infomap(&ig, NULL, NULL, num_trials, &membership, &codelength);
    }

    void clauset_newman_moore() {
        igraph_community_fastgreedy(&ig, NULL, &merges, &modularity, &membership);
    }

    void blondel_multilevel() {
        igraph_community_multilevel(&ig, NULL, &membership, NULL, NULL);
    }

    void walktrap_pons_latapy(int num_steps) {
        igraph_community_walktrap(&ig, NULL, num_steps, &merges, &modularity, &membership);
    }

    void label_propagation() {
        igraph_community_label_propagation(&ig, &membership, NULL, NULL, NULL, NULL);
    }

    void girvan_newman() {
        igraph_community_edge_betweenness(&ig, NULL, NULL, NULL, NULL, NULL, &membership, IGRAPH_UNDIRECTED, NULL);
    }

    int export_graph(CDLib::graph& g) {
        if (ig_status) {
            g.clear();
            for (int k = 0; k < igraph_vcount(&ig); k++)g.add_node();
            igraph_eit_t it;
            IGRAPH_CHECK(igraph_eit_create(&ig, igraph_ess_all(IGRAPH_EDGEORDER_FROM), &it));
            IGRAPH_FINALLY(igraph_eit_destroy, &it);
            while (!IGRAPH_EIT_END(it)) {
                igraph_integer_t from, to;
                igraph_edge(&ig, IGRAPH_EIT_GET(it), &from, &to);
                g.add_edge(static_cast<id_type> (from), static_cast<id_type> (to), 1);
                IGRAPH_EIT_NEXT(it);
            }
            g.set_graph_name(graph_name);
            igraph_eit_destroy(&it);
            IGRAPH_FINALLY_CLEAN(1);
            return 1;
        } else {
            return 0;
        }
    }

    ~igraph_ext() {
        igraph_destroy(&ig);
        igraph_vector_destroy(&membership);
        igraph_vector_destroy(&modularity);
        igraph_matrix_destroy(&merges);
    }
};

int igraph2graph(igraph_t& ig, CDLib::graph& g) {
    g.clear();
    for (int k = 0; k < igraph_vcount(&ig); k++)g.add_node();
    igraph_eit_t it;
    IGRAPH_CHECK(igraph_eit_create(&ig, igraph_ess_all(IGRAPH_EDGEORDER_FROM), &it));
    IGRAPH_FINALLY(igraph_eit_destroy, &it);
    while (!IGRAPH_EIT_END(it)) {
        igraph_integer_t from, to;
        igraph_edge(&ig, IGRAPH_EIT_GET(it), &from, &to);
        g.add_edge(static_cast<id_type> (from), static_cast<id_type> (to), 1);
        IGRAPH_EIT_NEXT(it);
    }
    igraph_eit_destroy(&it);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

int graph2igraph(const graph& g, igraph_t& ig) {
    igraph_vector_t edges = IGRAPH_VECTOR_NULL;
    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&edges, g.get_num_edges()));
    for (id_type i = 0; i < g.get_num_nodes(); i++) {
        for (adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i); aeit++) {
            if (i <= aeit->first) {
                IGRAPH_CHECK(igraph_vector_push_back(&edges, static_cast<long> (i)));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, static_cast<long> (aeit->first)));
            }
        }
    }
    IGRAPH_CHECK(igraph_create(&ig, &edges, g.get_num_nodes(), g.is_directed()));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return 1;
}

int graph2igraphw(const graph& g, igraph_t& ig, igraph_vector_t* wts) {
    igraph_vector_t edges = IGRAPH_VECTOR_NULL;
    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&edges, g.get_num_edges()));
    for (id_type i = 0; i < g.get_num_nodes(); i++) {
        for (adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i); aeit++) {
            if (i <= aeit->first) {
                IGRAPH_CHECK(igraph_vector_push_back(&edges, static_cast<long> (i)));
                IGRAPH_CHECK(igraph_vector_push_back(&edges, static_cast<long> (aeit->first)));
            }
        }
    }
    IGRAPH_CHECK(igraph_create(&ig, &edges, g.get_num_nodes(), g.is_directed()));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    igraph_vector_init(wts,igraph_ecount(&ig));
    for(igraph_integer_t i=0;i<igraph_ecount(&ig);i++){
        igraph_integer_t from,to;
        igraph_edge(&ig, i,&from, &to);
        VECTOR(*wts)[i] = g.get_edge_weight(from,to);
    }
    return 1;
}


#endif	/* IGRAPH_EXT_H */

