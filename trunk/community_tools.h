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
    
    
    double partition_quality(const graph& g,vector<node_set>& comms,double (*func)(const graph& g,node_set& comm));
    double modularity(const graph& g, vector<node_set>& comms);
    double modularity_density(const graph& g, vector<node_set>& comms);
    double community_score(const graph& g, vector<node_set>& comms);
    double description_length(const graph& g, vector<node_set>& comms);
    
    void compute_all_metrics_partition(const graph& g, vector<node_set>& comms,vector<double>& metrics);
    
    double rand_index(id_type num_nodes, vector<node_set>& comms1, vector<node_set>& comms2);
    double dongen_index(id_type num_nodes, vector<node_set>& comms1, vector<node_set>& comms2);
    double nmi(id_type num_nodes, vector<node_set>& comms1, vector<node_set>& comms2);
    double variation_of_information(id_type num_nodes, vector<node_set>& comms1, vector<node_set>& comms2);
    
    void convert_labels_to_communities(vector<id_type>& labels,vector<node_set>& communities);
    bool read_partition(const graph& g,const string& filepath,vector<node_set>& communities);
    
    void convert_communities_to_labels(const vector<node_set>& communities,vector<id_type>& labels);
    bool write_partition(const graph& g,const string& filepath,vector<node_set>& communities);
    bool write_partition_unlabelled(const graph& g,const string& filepath,vector<node_set>& communities);
    
    
    double evolutionary_cluster_validation(const graph& g1, const graph& g2,vector<node_set>& comms1,vector<node_set>& comms2,double (*func)(id_type num_nodes, vector<node_set>& comms1, vector<node_set>& comms2));
    double evolutionary_cluster_validation_union(const graph& g1, const graph& g2,vector<node_set>& comms1,vector<node_set>& comms2,double (*func)(id_type num_nodes, vector<node_set>& comms1, vector<node_set>& comms2));
    
    bool is_member_of(id_type id, node_set& ns);
    bool in_same_comm(id_type i,id_type j,vector<node_set>& comms);
    id_type in_comm(id_type i, vector<node_set>& comms);
    double degree_homogenity_test(graph& g, node_set & ns);
    double kl_divergence(vector<double>& p,vector<double>& q);
    double entropy_comparision_test(const graph& g,node_set& ns);
    
    void get_community_graph(const graph&g, vector<node_set>& comms,graph& comm_graph);
    
    void compute_confusion_matrix_local(const graph& g,node_set& observed, node_set& truth,confusion_matrix_local& res);
};

#endif	/* COMMUNITY_TOOLS_H */

