/*
 * File:   community_tools.h
 * Author: bharath
 *
 * Created on 5 April, 2012, 2:27 AM
 */

#ifndef COMMUNITY_TOOLS_H
#define	COMMUNITY_TOOLS_H

#include "graph.h"
#include "graph_operations.h"
#include "paths_and_components.h"
#include "statistics.h"

using namespace std;

namespace CDLib
{

    struct confusion_matrix_local
    {
        id_type true_positives;
        id_type false_positives;
        id_type true_negatives;
        id_type false_negatives;
        double TPR;
        double FPR;
        double FNR;
        double TNR;
        double wTPR;
        double wFPR;
        double precision;
        double recall;
        double f_score;
    };

    struct community_metrics
    {
        id_type size;
        double intercluster_edges;
        double intracluster_edges;
        double internal_density;
        double conductance;
        double nassoc;
        double cohesion;
        double expansion;
        double modularity;
        double modularity_density;
        double community_score;
        double degree_homogenity;
        double degree_entropy;
        double rwalk_entropy;
    };

    struct cluster_edges
    {
        id_type num_inter_cluster_edges;
        id_type num_intra_cluster_edges;

        double num_expected_inter_cluster_edges;
        double num_expected_intra_cluster_edges;

        double wt_inter_cluster_edges;
        double wt_intra_cluster_edges;

        double wt_expected_inter_cluster_edges;
        double wt_expected_intra_cluster_edges;

        double max_odf;
        double avg_odf;
        double flake_odf;
        bool is_weak_radicchi_community;
        bool is_strong_radicchi_community;
        cluster_edges();
        cluster_edges(const graph& g, node_set& comm);
    };

    double volume_comm(const graph& g, node_set& comm);
    double cut_comm(const graph& g, node_set& comm);
    double ratio_cut_comm(const graph& g, node_set& comm);
    double conductance_comm(const graph& g, node_set& comm);
    double resistance_comm(const graph& g, node_set& comm);
    double normalized_cut_comm(const graph& g, node_set& comm);
    double expansion_comm(const graph& g, node_set& comm);
    double internal_density_comm(const graph& g, node_set& comm);
    double max_odf_comm(const graph& g, node_set& comm);
    double avg_odf_comm(const graph& g, node_set& comm);
    double flake_odf_comm(const graph& g, node_set& comm);
    bool radicchi_community(const graph& g, node_set& comm,bool strong);
    double modularity_comm(const graph& g, node_set& comm);
    double modularity_density_comm(const graph& g, node_set& comm);
    double ratio_assoc_comm(const graph& g, node_set& comm);
    double normalized_assoc_comm(const graph& g, node_set& comm);

    double partition_quality(const graph& g,vector<node_set>& comms,double (*func)(const graph& g,node_set& comm));
    double modularity(const graph& g, vector<node_set>& comms);
    double modularity_density(const graph& g, vector<node_set>& comms);
    double community_score(const graph& g, vector<node_set>& comms);
    double description_length(const graph& g, vector<node_set>& comms);
    void compute_community_metrics(const graph& g, node_set& comm,community_metrics& metrics);
    void compute_all_metrics_partition(const graph& g, vector<node_set>& comms,vector<double>& metrics);
    void compute_imp_metrics_partition(const graph& g, vector<node_set>& comms,vector<double>& metrics);
    double rand_index(id_type num_nodes, vector<node_set>& comms1, vector<node_set>& comms2);
    double dongen_index(id_type num_nodes, vector<node_set>& comms1, vector<node_set>& comms2);
    double nmi(id_type num_nodes, vector<node_set>& comms1, vector<node_set>& comms2);
    double variation_of_information(id_type num_nodes, vector<node_set>& comms1, vector<node_set>& comms2);
    id_type get_num_of_communities(const vector<id_type>& labels);
    void convert_labels_to_communities(const vector<id_type>& labels,vector<node_set>& communities);
    bool read_partition(const graph& g,const string& filepath,vector<node_set>& communities);

    void convert_communities_to_labels(const  vector<node_set>& communities,vector<id_type>& labels);

    id_type reindex_communities(const vector<id_type>& old_comms,vector<id_type>& new_comms);
    id_type reindex_communities(vector<id_type>& labels);

    bool write_partition(const graph& g,const string& filepath,vector<node_set>& communities);
    bool write_partition_unlabelled(const graph& g,const string& filepath,vector<node_set>& communities);


    double evolutionary_cluster_validation(const graph& g1, const graph& g2,vector<node_set>& comms1,vector<node_set>& comms2,double (*func)(id_type num_nodes, vector<node_set>& comms1, vector<node_set>& comms2));
    double evolutionary_cluster_validation_union(const graph& g1, const graph& g2,vector<node_set>& comms1,vector<node_set>& comms2,double (*func)(id_type num_nodes, vector<node_set>& comms1, vector<node_set>& comms2));

    bool is_member_of(id_type id, node_set& ns);
    bool in_same_comm(id_type i,id_type j,vector<node_set>& comms);
    id_type in_comm(id_type i, vector<node_set>& comms);
    double degree_homogenity_test(graph& g, node_set & ns);
    double entropy_comparision_test(const graph& g,node_set& ns);

    void get_community_graph(const graph&g, vector<node_set>& comms,graph& comm_graph);

    void compute_confusion_matrix_local(const graph& g,node_set& observed, node_set& truth,confusion_matrix_local& res);

    void compute_confusion_matrix_global(const graph&g, vector<node_set>& observed,vector<node_set>& truth,vector< vector<id_type> >& cmat);
    void componentize_and_reindex_labels(const graph& g,const vector<id_type>& templabels, vector<id_type>& labels);
    void componentize_and_reindex_labels(const graph& g,vector<id_type>& labels);
};

#endif	/* COMMUNITY_TOOLS_H */

