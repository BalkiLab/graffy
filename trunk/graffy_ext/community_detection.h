/* 
 * File:   community_detection.h
 * Author: bharath
 *
 * Created on 24 April, 2013, 6:17 AM
 */

#ifndef COMMUNITY_DETECTION_H
#define	COMMUNITY_DETECTION_H

#include "includes.h"
#include "igraph.h"
#include "xmeans_ext.h"
#include "igraph_ext.h"

namespace CDA {

    void louvain_method_igraph_weighted(const graph& g, double t, const vector<id_type>& init_partition, vector<id_type>& labels) {
        igraph_t ig;
        igraph_vector_t ilabels, flabels, wts;
        graph2igraphw(g, ig, &wts);
        igraph_vector_init(&flabels, 0);
        igraph_vector_init(&ilabels, g.get_num_nodes());
        if (init_partition.size() == g.get_num_nodes()) {
            for (id_type i = 0; i < init_partition.size(); i++)
                VECTOR(ilabels)[i] = init_partition[i];
            igraph_community_multilevel_stability_incremental(&ig, &wts, t, &ilabels, &flabels, NULL, NULL);
        } else igraph_community_multilevel_stability(&ig, &wts, t, &flabels, NULL, NULL);
        labels.assign(g.get_num_nodes(), 0);
        for (id_type i = 0; i < g.get_num_nodes(); i++) labels[i] = VECTOR(flabels)[i];
        igraph_vector_destroy(&flabels);
        igraph_vector_destroy(&ilabels);
        igraph_vector_destroy(&wts);
        igraph_destroy(&ig);
        CDLib::componentize_and_reindex_labels(g, labels);
    }

    inline void louvain_method_igraph_weighted(const graph& g, double t, vector<id_type>& labels) {
        louvain_method_igraph_weighted(g, t, vector<id_type>(), labels);
    }

    void wlogv_igraph_weighted(const graph& g, vector<id_type>& labels) {
        igraph_t ig;
        igraph_vector_t flabels, wts;
        graph2igraphw(g, ig, &wts);
        igraph_vector_init(&flabels, 0);
        igraph_community_multilevel_wlogv(&ig, &wts, &flabels, NULL, NULL);
        labels.assign(g.get_num_nodes(), 0);
        for (id_type i = 0; i < g.get_num_nodes(); i++) labels[i] = VECTOR(flabels)[i];
        igraph_vector_destroy(&flabels);
        igraph_vector_destroy(&wts);
        igraph_destroy(&ig);
        CDLib::componentize_and_reindex_labels(g, labels);
    }

    void louvain_method_igraph(const graph& g, double t, const vector<id_type>& init_partition, vector<id_type>& labels) {
        igraph_t ig;
        graph2igraph(g, ig);
        igraph_vector_t ilabels, flabels;
        igraph_vector_init(&flabels, 0);
        igraph_vector_init(&ilabels, g.get_num_nodes());
        if (init_partition.size() == g.get_num_nodes()) {
            for (id_type i = 0; i < init_partition.size(); i++)
                VECTOR(ilabels)[i] = init_partition[i];
            igraph_community_multilevel_stability_incremental(&ig, NULL, t, &ilabels, &flabels, NULL, NULL);
        } else igraph_community_multilevel_stability(&ig, NULL, t, &flabels, NULL, NULL);
        labels.assign(g.get_num_nodes(), 0);
        for (id_type i = 0; i < g.get_num_nodes(); i++) labels[i] = VECTOR(flabels)[i];
        igraph_vector_destroy(&flabels);
        igraph_vector_destroy(&ilabels);
        igraph_destroy(&ig);
        CDLib::componentize_and_reindex_labels(g, labels);
    }

    inline void louvain_method_igraph(const graph& g, double t, vector<id_type>& labels) {
        louvain_method_igraph(g, t, vector<id_type>(), labels);

    }

    void label_propagation_igraph(const graph& g, const vector<id_type>& init_partition, vector<id_type>& labels) {
        igraph_t ig;
        graph2igraph(g, ig);
        igraph_vector_t ilabels, flabels;
        igraph_vector_init(&flabels, 0);
        igraph_vector_init(&ilabels, g.get_num_nodes());
        if (init_partition.size() == g.get_num_nodes()) {
            for (id_type i = 0; i < init_partition.size(); i++)
                VECTOR(ilabels)[i] = init_partition[i];
            igraph_community_label_propagation(&ig, &flabels, NULL, &ilabels, NULL, NULL);
        } else igraph_community_label_propagation(&ig, &flabels, NULL, NULL, NULL, NULL);
        labels.assign(g.get_num_nodes(), 0);
        for (id_type i = 0; i < g.get_num_nodes(); i++) labels[i] = VECTOR(flabels)[i];
        igraph_vector_destroy(&flabels);
        igraph_vector_destroy(&ilabels);
        igraph_destroy(&ig);
        CDLib::componentize_and_reindex_labels(g, labels);
    }

    inline void label_propagation_igraph(const graph&g, vector<id_type>& labels) {
        label_propagation_igraph(g, vector<id_type>(), labels);
    }

    void infomap_igraph(const graph& g, vector<id_type>& labels) {
        igraph_t ig;
        igraph_real_t cl;
        graph2igraph(g, ig);
        igraph_vector_t flabels;
        igraph_vector_init(&flabels, 0);
        igraph_community_infomap(&ig, NULL, NULL, 1, &flabels, &cl);
        labels.assign(g.get_num_nodes(), 0);
        for (id_type i = 0; i < g.get_num_nodes(); i++) labels[i] = VECTOR(flabels)[i];
        igraph_vector_destroy(&flabels);
        CDLib::componentize_and_reindex_labels(g, labels);
    }

    void wlogv_igraph(const graph& g, vector<id_type>& labels) {
        igraph_t ig;
        graph2igraph(g, ig);
        igraph_vector_t flabels;
        igraph_vector_init(&flabels, 0);
        igraph_community_multilevel_wlogv(&ig, NULL, &flabels, NULL, NULL);
        labels.assign(g.get_num_nodes(), 0);
        for (id_type i = 0; i < g.get_num_nodes(); i++) labels[i] = VECTOR(flabels)[i];
        igraph_vector_destroy(&flabels);
        CDLib::componentize_and_reindex_labels(g, labels);
    }

    void wlogv2_igraph(const graph& g, vector<id_type>& labels) {
        igraph_t ig;
        graph2igraph(g, ig);
        igraph_vector_t flabels;
        igraph_vector_init(&flabels, 0);
        igraph_community_multilevel_wlogv2(&ig, NULL, &flabels, NULL, NULL);
        labels.assign(g.get_num_nodes(), 0);
        for (id_type i = 0; i < g.get_num_nodes(); i++) labels[i] = VECTOR(flabels)[i];
        igraph_vector_destroy(&flabels);
        CDLib::componentize_and_reindex_labels(g, labels);
    }

    void parabola_igraph(const graph& g, vector<id_type>& labels) {
        igraph_t ig;
        graph2igraph(g, ig);
        igraph_vector_t flabels;
        igraph_vector_init(&flabels, 0);
        igraph_community_multilevel_parabola(&ig, NULL, &flabels, NULL, NULL);
        labels.assign(g.get_num_nodes(), 0);
        for (id_type i = 0; i < g.get_num_nodes(); i++) labels[i] = VECTOR(flabels)[i];
        igraph_vector_destroy(&flabels);
        CDLib::componentize_and_reindex_labels(g, labels);
    }

    void vertex_mover_igraph(const graph& g, double t, const vector<id_type>& init_partition, vector<id_type>& labels) {
        igraph_t ig;
        graph2igraph(g, ig);
        igraph_vector_t ilabels, flabels;
        igraph_vector_init(&flabels, 0);
        igraph_vector_init(&ilabels, g.get_num_nodes());
        if (init_partition.size() == g.get_num_nodes()) {
            for (id_type i = 0; i < init_partition.size(); i++)
                VECTOR(ilabels)[i] = init_partition[i];
            igraph_i_community_vertex_mover(&ig, NULL, &ilabels, t, &flabels, NULL);
        } else igraph_i_community_vertex_mover(&ig, NULL, NULL, t, &flabels, NULL);
        labels.assign(g.get_num_nodes(), 0);
        for (id_type i = 0; i < g.get_num_nodes(); i++) labels[i] = VECTOR(flabels)[i];
        igraph_vector_destroy(&flabels);
        igraph_vector_destroy(&ilabels);
        igraph_destroy(&ig);
        CDLib::componentize_and_reindex_labels(g, labels);
    }

    inline void vertex_mover_igraph(const graph&g, double t, vector<id_type>& labels) {
        vertex_mover_igraph(g, t, vector<id_type>(), labels);
    }

    void vertex_mover_igraph_wlogv(const graph& g, const vector<id_type>& init_partition, vector<id_type>& labels) {
        igraph_t ig;
        graph2igraph(g, ig);
        igraph_vector_t ilabels, flabels;
        igraph_vector_init(&flabels, 0);
        igraph_vector_init(&ilabels, g.get_num_nodes());
        if (init_partition.size() == g.get_num_nodes()) {
            for (id_type i = 0; i < init_partition.size(); i++)
                VECTOR(ilabels)[i] = init_partition[i];
            igraph_i_community_vertex_mover_wlogv(&ig, NULL, &ilabels, &flabels, NULL);
        } else igraph_i_community_vertex_mover_wlogv(&ig, NULL, NULL, &flabels, NULL);
        labels.assign(g.get_num_nodes(), 0);
        for (id_type i = 0; i < g.get_num_nodes(); i++) labels[i] = VECTOR(flabels)[i];
        igraph_vector_destroy(&flabels);
        igraph_vector_destroy(&ilabels);
        igraph_destroy(&ig);
        CDLib::componentize_and_reindex_labels(g, labels);
    }

    inline void vertex_mover_igraph_wlogv(const graph&g, vector<id_type>& labels) {
        vertex_mover_igraph_wlogv(g, vector<id_type>(), labels);
    }

    void vertex_mover_igraph_weighted(const graph& g, double t, const vector<id_type>& init_partition, vector<id_type>& labels) {
        igraph_t ig;
        igraph_vector_t ilabels, flabels, wts;
        graph2igraphw(g, ig, &wts);
        igraph_vector_init(&flabels, 0);
        igraph_vector_init(&ilabels, g.get_num_nodes());
        if (init_partition.size() == g.get_num_nodes()) {
            for (id_type i = 0; i < init_partition.size(); i++)
                VECTOR(ilabels)[i] = init_partition[i];
            igraph_i_community_vertex_mover(&ig, &wts, &ilabels, t, &flabels, NULL);
        } else igraph_i_community_vertex_mover(&ig, &wts, NULL, t, &flabels, NULL);
        labels.assign(g.get_num_nodes(), 0);
        for (id_type i = 0; i < g.get_num_nodes(); i++) labels[i] = VECTOR(flabels)[i];
        igraph_vector_destroy(&flabels);
        igraph_vector_destroy(&ilabels);
        igraph_vector_destroy(&wts);
        igraph_destroy(&ig);
        CDLib::componentize_and_reindex_labels(g, labels);
    }

    void vertex_mover_igraph_weighted(const graph& g, double t, vector<id_type>& labels) {
        vertex_mover_igraph_weighted(g, t, vector<id_type>(), labels);
    }

    void vertex_mover_igraph_weighted_wlogv(const graph& g, const vector<id_type>& init_partition, vector<id_type>& labels) {
        igraph_t ig;
        igraph_vector_t ilabels, flabels, wts;
        graph2igraphw(g, ig, &wts);
        igraph_vector_init(&flabels, 0);
        igraph_vector_init(&ilabels, g.get_num_nodes());
        if (init_partition.size() == g.get_num_nodes()) {
            for (id_type i = 0; i < init_partition.size(); i++)
                VECTOR(ilabels)[i] = init_partition[i];
            igraph_i_community_vertex_mover_wlogv(&ig, &wts, &ilabels, &flabels, NULL);
        } else igraph_i_community_vertex_mover_wlogv(&ig, &wts, NULL, &flabels, NULL);
        labels.assign(g.get_num_nodes(), 0);
        for (id_type i = 0; i < g.get_num_nodes(); i++) labels[i] = VECTOR(flabels)[i];
        igraph_vector_destroy(&flabels);
        igraph_vector_destroy(&ilabels);
        igraph_vector_destroy(&wts);
        igraph_destroy(&ig);
        CDLib::componentize_and_reindex_labels(g, labels);
    }

    void vertex_mover_igraph_weighted_wlogv(const graph& g, vector<id_type>& labels) {
        vertex_mover_igraph_weighted_wlogv(g, vector<id_type>(), labels);
    }

    void igraph_walktrap(const graph& g, int num_steps, vector<id_type>& labels) {

        igraph_t ig;
        graph2igraph(g, ig);
        igraph_vector_t flabels;
        igraph_vector_init(&flabels, 0);
        igraph_community_walktrap(&ig, NULL, num_steps, NULL, NULL, &flabels);
        labels.assign(g.get_num_nodes(), 0);
        for (id_type i = 0; i < g.get_num_nodes(); i++) labels[i] = VECTOR(flabels)[i];
        igraph_vector_destroy(&flabels);
        CDLib::componentize_and_reindex_labels(g, labels);

    }

    void k_core_igraph(const graph& g, id_type k, node_set& core_nodes) {
        igraph_t ig;
        graph2igraph(g, ig);
        igraph_vector_t cores;
        igraph_vector_init(&cores, 0);
        igraph_coreness(&ig, &cores, IGRAPH_ALL);
        for (id_type i = 0; i < static_cast<id_type> (igraph_vector_size(&cores)); i++) {
            if (VECTOR(cores)[i] >= k) core_nodes.insert(i);
        }
        igraph_vector_destroy(&cores);
        igraph_destroy(&ig);
    }

    int mqi_single(const graph& g, const node_set& ns_in, node_set& comms_out) {
        id_type ctr = 3;
        double a = 0, c = 0;
        unordered_map<string, id_type> newnodes;
        vector<string> newnodes_rm(3, "");
        deque<pair<id_type, pair<id_type, double> > > newedges;
        for (node_set::const_iterator nit = ns_in.begin(); nit != ns_in.end(); nit++) {
            newnodes.insert(make_pair(g.get_node_label(*nit), ctr++));
            newnodes_rm.push_back(g.get_node_label(*nit));
            for (adjacent_edges_iterator aeit = g.out_edges_begin(*nit); aeit != g.out_edges_end(*nit); aeit++)
                if (ns_in.find(aeit->first) == ns_in.end())
                    c += aeit->second;
        }
        a = ns_in.size();
        for (node_set::const_iterator nit = ns_in.begin(); nit != ns_in.end(); nit++) {
            string to_node_label = g.get_node_label(*nit);
            id_type to_node_id = newnodes[to_node_label];
            for (adjacent_edges_iterator aeit = g.in_edges_begin(*nit); aeit != g.in_edges_end(*nit); aeit++) {
                if ((ns_in.find(aeit->first) != ns_in.end()) && (*nit != aeit->first)) {
                    string from_node_label = g.get_node_label(aeit->first);
                    id_type from_node_id = newnodes[from_node_label];
                    newedges.push_back(make_pair(to_node_id, make_pair(from_node_id, a)));
                } else newedges.push_back(make_pair(1, make_pair(to_node_id, a)));
            }
            newedges.push_back(make_pair(to_node_id, make_pair(2, c)));
        }
        string outfilename(tmpnam(NULL)), resfilename(tmpnam(NULL));
        ofstream outfile(outfilename);
        outfile << "p max " << newnodes.size() + 2 << " " << newedges.size() << endl;
        outfile << "n 1 s" << endl << "n 2 t" << endl;
        for (id_type k = 0; k < newedges.size(); k++) outfile << "a " << newedges[k].first << " " << newedges[k].second.first << " " << newedges[k].second.second << endl;
        outfile.close();
        ostringstream oss, oss2;
        oss << "hi_pr <" << outfilename + " | grep -e '^c [0-9].*' | sed -e 's/c //g' > " << resfilename;
        int r1 = system(oss.str().c_str());
        ifstream ifs(resfilename);
        node_set wrong;
        while (!ifs.eof()) {
            id_type temp;
            ifs >> temp;
            if (temp > 2) wrong.insert(g.get_node_id(newnodes_rm[temp]));
        }
        ifs.close();
        oss2 << "rm " << outfilename + " " << resfilename;
        int r2 = system(oss2.str().c_str());
        for (id_type i = 3; i < newnodes_rm.size(); i++) {
            id_type act_node_id = g.get_node_id(newnodes_rm[i]);
            if (!is_member_of(act_node_id, wrong)) comms_out.insert(act_node_id);
        }
        return r1 && r2;
    }

    id_type mqi_feedback(const graph& g, const node_set& comm, id_type num_iters, node_set& sideA, node_set& sideB) {
        node_set last_comm = comm;
        id_type i = 0, last_comm_size = 0;
        for (; i < num_iters; i++) {
            last_comm_size = last_comm.size();
            mqi_single(g, last_comm, sideA);
            if (sideA.size() == last_comm_size) break;
            last_comm = sideA;
        }
        for (node_set::const_iterator nit = comm.begin(); nit != comm.end(); nit++)
            if (sideA.find(*nit) == sideA.end())sideB.insert(*nit);
        return i;
    }

    void botgrep_prefiltering(const graph &g, const vector<double>& q, double r, id_type max_k, vector<id_type>& final_labels) {
        vector<double> qout;
        timer_rt op_timer;
        op_timer.start_clock();
        run_random_walks(g, q, log2(g.get_num_nodes()), 0, qout);
        op_timer.stop_clock();
        cerr << "Completed Random Walks in " << op_timer.run_time() << endl;
        vector< vector<double> > s(g.get_num_nodes(), vector<double>(1, 0));
        op_timer.start_clock();
        for (id_type i = 0; i < g.get_num_nodes(); i++) {
            if (g.get_node_out_weight(i)) {
                s[i][0] = pow(qout[i] / g.get_node_out_weight(i), 1 / r);
            }
        }
        vector< vector<double> > centroids;
        XMeans::xmeans(s, XMeans::xmeans_options(2, max_k), final_labels, centroids);
        op_timer.stop_clock();
        cerr << "Completed Xmeans For r=" << r << " maxk=" << max_k << " in " << op_timer.run_time() << endl;
        CDLib::componentize_and_reindex_labels(g, final_labels);
    }

    void botgrep_prefiltering(const graph &g, double r, id_type max_k, vector<id_type>& final_labels) {
        vector<double> q(g.get_num_nodes(), 1.0 / g.get_num_nodes());
        botgrep_prefiltering(g, q, r, max_k, final_labels);
    }

    int botgrep_sybillinfer(const graph& g, const node_set true_bots_in_comm, const node_set& non_bots_in_comm,id_type num_total_bot,id_type num_total_benign, double conductance, id_type max_iters, const string& outfile) {
        node_set bots, not_bots;
        bots.clear();
        not_bots.clear();
        graph gen(0, 0);
        vector<id_type> ids;
        for (node_set::const_iterator nit = true_bots_in_comm.begin(); nit != true_bots_in_comm.end(); nit++) {
            ids.push_back(*nit);
            gen.add_node();
        }
        for (node_set::const_iterator nit = non_bots_in_comm.begin(); nit != non_bots_in_comm.end(); nit++) {
            ids.push_back(*nit);
            gen.add_node();
        }
        unordered_map<id_type, id_type> rmap;
        for (id_type i = 0; i < ids.size(); i++)rmap.insert(make_pair(ids[i], i));
        for (unordered_map<id_type, id_type>::iterator mapit = rmap.begin(); mapit != rmap.end(); mapit++) {
            id_type node_id = mapit->second;
            for (adjacent_edges_iterator aeit = g.out_edges_begin(mapit->first); aeit != g.out_edges_end(mapit->first); aeit++) {
                unordered_map<id_type, id_type>::iterator ret = rmap.find(aeit->first);
                if (rmap.find(aeit->first) != rmap.end()) {
                    id_type neigh_id = ret->second;
                    gen.add_edge(node_id, neigh_id, 1);
                }
            }
        }
        string tempfile = outfile;
        ofstream ofs(tempfile + ".si_input");
        for (id_type i = 0; i < gen.get_num_nodes(); i++) {
            deque<id_type>edges;
            for (adjacent_edges_iterator aeit = gen.out_edges_begin(i); aeit != gen.out_edges_end(i); aeit++) edges.push_back(aeit->first);
            sort(edges.begin(), edges.end());
            ofs << i << " " << ids[i] << " ";
            for (id_type j = 0; j < edges.size(); j++) ofs << edges[j] << " ";
            ofs << endl;
        }
        ofs.close();
        ostringstream oss;
        oss << "/opt/Extensions/modified_sybillinfer/botinfer_multi_round " << gen.get_num_nodes() << " " << true_bots_in_comm.size() << " " << num_total_benign<< " " << num_total_bot<< " " << time(NULL) << " " << tempfile << ".si_input " << conductance << " " << max_iters << " " << tempfile << " 2>>errors.log >" << outfile << "-si-res.txt";
        int r1 = system(oss.str().c_str());
        ifstream ifs(tempfile + ".list");
        while (!ifs.eof()) {
            id_type temp;
            ifs >> temp;
            unordered_map<id_type, id_type>::const_iterator ret = rmap.find(temp);
            if (ret != rmap.end())bots.insert(temp);
        }
        for (unordered_map<id_type, id_type>::const_iterator it = rmap.begin(); it != rmap.end(); it++) {
            if (bots.find(it->first) == bots.end()) not_bots.insert(it->first);
        }
        ifs.close();
        string cmd = "rm " + tempfile + ".list " + tempfile + ".si_input";
        int r2 = system(cmd.c_str());
        return r1 && r2;
    }

    void louvain_igraph_generic(const graph&g, const string& algo, double t, vector<id_type>& vmlabels, vector<id_type>& labels) {
        igraph_t ig;
        igraph_vector_t flabels, wts;
        igraph_matrix_t memberships;
        igraph_matrix_init(&memberships, 0, 0);
        igraph_vector_t* wt_ptr = NULL;
        bool is_weighted = g.is_weighted();
        if (is_weighted) {
            graph2igraphw(g, ig, &wts);
            wt_ptr = &wts;
        } else graph2igraph(g, ig);
        igraph_vector_init(&flabels, 0);
        if (algo == "wlogv") igraph_community_multilevel_wlogv(&ig, wt_ptr, &flabels, &memberships, NULL);
        else igraph_community_multilevel_stability(&ig, wt_ptr, t, &flabels, &memberships, NULL);
        labels.assign(g.get_num_nodes(), 0);
        vmlabels.assign(g.get_num_nodes(), 0);
        for (id_type i = 0; i < g.get_num_nodes(); i++) {
            labels[i] = VECTOR(flabels)[i];
            vmlabels[i] = MATRIX(memberships, 0, i);
        }
        igraph_vector_destroy(&flabels);
        igraph_matrix_destroy(&memberships);
        if (g.is_weighted()) igraph_vector_destroy(&wts);
        igraph_destroy(&ig);
        CDLib::componentize_and_reindex_labels(g, labels);
    }

    void leading_eigenvector_igraph(const graph& g, vector<id_type>& labels) {
        igraph_t ig;
        igraph_arpack_options_t options;
        igraph_arpack_options_init(&options);
        igraph_vector_t flabels;
        graph2igraph(g, ig);
        igraph_vector_init(&flabels, 0);
        igraph_community_leading_eigenvector(&ig, NULL, &flabels, g.get_num_nodes(), &options, NULL, 0, NULL, NULL,NULL,NULL,NULL);
        labels.assign(g.get_num_nodes(), 0);
        for (id_type i = 0; i < g.get_num_nodes(); i++) labels[i] = VECTOR(flabels)[i];
        igraph_vector_destroy(&flabels);
        igraph_destroy(&ig);
        CDLib::componentize_and_reindex_labels(g, labels);
    }

    void cnm_igraph_weighted(const graph& g, vector<id_type>& labels) {
        igraph_t ig;
        igraph_vector_t flabels, wts;
        graph2igraphw(g, ig, &wts);
        igraph_vector_init(&flabels, 0);
        igraph_community_fastgreedy(&ig, &wts, NULL, NULL, &flabels);
        labels.assign(g.get_num_nodes(), 0);
        for (id_type i = 0; i < g.get_num_nodes(); i++) labels[i] = VECTOR(flabels)[i];
        igraph_vector_destroy(&flabels);
        igraph_vector_destroy(&wts);
        igraph_destroy(&ig);
        CDLib::componentize_and_reindex_labels(g, labels);
    }

    void cnm_igraph(const graph& g, vector<id_type>& labels) {
        igraph_t ig;
        igraph_vector_t flabels;
        graph2igraph(g, ig);
        igraph_vector_init(&flabels, 0);
        igraph_community_fastgreedy(&ig, NULL, NULL, NULL, &flabels);
        labels.assign(g.get_num_nodes(), 0);
        for (id_type i = 0; i < g.get_num_nodes(); i++) labels[i] = VECTOR(flabels)[i];
        igraph_vector_destroy(&flabels);
        igraph_destroy(&ig);
        CDLib::componentize_and_reindex_labels(g, labels);
    }
    
        void gn_igraph_weighted(const graph& g, vector<id_type>& labels) {
        igraph_t ig;
        igraph_vector_t flabels, wts;
        graph2igraphw(g, ig, &wts);
        igraph_vector_init(&flabels, 0);
        igraph_community_edge_betweenness(&ig,NULL,NULL,NULL,NULL,NULL,&flabels,g.is_directed(),&wts);
        labels.assign(g.get_num_nodes(), 0);
        for (id_type i = 0; i < g.get_num_nodes(); i++) labels[i] = VECTOR(flabels)[i];
        igraph_vector_destroy(&flabels);
        igraph_vector_destroy(&wts);
        igraph_destroy(&ig);
        CDLib::componentize_and_reindex_labels(g, labels);
    }

    void gn_igraph(const graph& g, vector<id_type>& labels) {
        igraph_t ig;
        igraph_vector_t flabels;
        graph2igraph(g, ig);
        igraph_vector_init(&flabels, 0);
        igraph_community_edge_betweenness(&ig,NULL,NULL,NULL,NULL,NULL,&flabels,g.is_directed(),NULL);
        labels.assign(g.get_num_nodes(), 0);
        for (id_type i = 0; i < g.get_num_nodes(); i++) labels[i] = VECTOR(flabels)[i];
        igraph_vector_destroy(&flabels);
        igraph_destroy(&ig);
        CDLib::componentize_and_reindex_labels(g, labels);
    }
    


};

#endif	/* COMMUNITY_DETECTION_H */

