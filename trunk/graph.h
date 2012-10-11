/* 
 * File:   graph.h
 * Author: bharath
 *
 * Created on April 1, 2012, 12:00 PM
 */

#ifndef GRAPH_H
#define	GRAPH_H

#include "bidirectional_label_map.h"
#include "double_adjacency_map.h"
#include "typedefs.h"
#include "utility.h"

namespace CDLib
{
    class graph {
    private:
        
        bool b_directed;
        bool b_weighted;
        bidirectional_label_map blm_labels;
        double_adjacency_map dam_backend;
        
    public:
        
        graph();
        graph(bool directed, bool weighted);
        
        bool is_directed() const;
        bool is_weighted() const;
        
        id_type get_num_nodes() const;
        id_type get_num_edges() const;
        id_type get_num_self_edges() const;
        wt_t get_total_weight() const;
        wt_t get_self_edges_weight() const;
        wt_t get_density() const;
         
        string get_node_label(id_type id) const;
        id_type get_node_id(const string& label) const;
        
        id_type get_node_in_degree(id_type id) const;
        id_type get_node_in_degree(const string& label) const; 
        
        id_type get_node_out_degree(id_type id) const;
        id_type get_node_out_degree(const string& label) const; 
        
        wt_t get_node_in_weight(id_type) const;
        wt_t get_node_in_weight(const string& label) const;
        
        wt_t get_node_out_weight(id_type) const;
        wt_t get_node_out_weight(const string& label) const;
        
        wt_t get_edge_weight(id_type from_id, id_type to_id) const;
        wt_t get_edge_weight(const string& from_label, const string& to_label) const;
        
        node_label_iterator node_labels_begin() const;
        node_label_iterator node_labels_end() const;
        
        adjacent_edges_iterator in_edges_begin(id_type id) const;
        adjacent_edges_iterator in_edges_begin(const string& label) const;

        adjacent_edges_iterator in_edges_end(id_type id) const;
        adjacent_edges_iterator in_edges_end(const string& label) const;

        adjacent_edges_iterator out_edges_begin(id_type id) const;
        adjacent_edges_iterator out_edges_begin(const string& label) const;

        adjacent_edges_iterator out_edges_end(id_type id) const;
        adjacent_edges_iterator out_edges_end(const string& label) const;
        
        id_type add_node(const string& label);
        id_type add_node();
        
        bool add_edge(id_type from_id, id_type to_id, wt_t weight);
        void add_self_edges(double weight);
        void remove_self_edges();
        wt_t add_edge(const string& from_label, const string& to_label,wt_t weight);
        
        bool remove_node(id_type id);
        bool remove_node(const string& label);
        id_type remove_isolates();
        
        bool remove_edge(id_type from_id, id_type to_id);
        wt_t remove_edge(const string& from_label, const string& to_label);
        
        bool set_edge_weight(id_type from_id, id_type to_id, wt_t weight);
        wt_t set_edge_weight(const string& from_label, const string& to_label,wt_t weight);
        
        bool remove_all_edges();
        bool clear();
        
        bool convert_to_unweighted(double threshold);
        bool convert_to_undirected();
        bool convert_to_directed();
        bool convert_to_weighted();
        
        double extreme_weight(bool max) const;
        double minimum_weight() const;
        double maximum_weight() const;
    };
};

#endif	/* GRAPH_H */

