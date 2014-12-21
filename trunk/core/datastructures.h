/* 
 * File:   datastructures.h
 * Author: bharath
 *
 * Created on 31 January, 2013, 8:16 PM
 */

#ifndef DATASTRUCTURES_H
#define	DATASTRUCTURES_H

#include "typedefs.h"
using namespace std;

namespace CDLib {

    typedef unordered_map<string, id_type> forward_map;
    typedef unordered_map<id_type, string> reverse_map;
    typedef forward_map::const_iterator node_label_iterator;

    class bidirectional_label_map {
    private:
        forward_map fm_labels;
        reverse_map rm_ids;
    public:
        bool swap_labels(id_type old_id,id_type new_id);
        id_type size() const;
        node_label_iterator begin() const;
        node_label_iterator end() const;
        string get_label(id_type id) const;
        id_type get_id(const string& label) const;
        bool insert(const string& label);
        bool erase(const string& label);
        bool erase(id_type id);
        bool clear();
    };

    typedef unordered_map<id_type, wt_t> adjacent_edge_sequence;
    typedef vector<adjacent_edge_sequence> adjacency_map;
    typedef typename adjacent_edge_sequence::const_iterator adjacent_edges_iterator;

    class double_adjacency_map {
    private:
        id_type st_num_edges;
        id_type st_num_self_edges;
        wt_t wt_total_wt;
        wt_t wt_self_edge_wt;
        adjacency_map am_in_edges;
        adjacency_map am_out_edges;
        vector<wt_t> vw_in_degree;
        vector<wt_t> vw_out_degree;
    public:
        double_adjacency_map();
        id_type num_nodes() const;
        id_type num_edges() const;
        id_type num_self_edges() const;
        wt_t total_weight() const;
        wt_t self_edges_weight() const;
        bool is_valid_node(id_type id) const;
        id_type in_degree(id_type id) const;
        id_type out_degree(id_type id) const;
        wt_t in_degree_wt(id_type id) const;
        wt_t out_degree_wt(id_type id) const;
        adjacent_edges_iterator in_edges_begin(id_type id) const;
        adjacent_edges_iterator in_edges_end(id_type id) const;
        adjacent_edges_iterator out_edges_begin(id_type id) const;
        adjacent_edges_iterator out_edges_end(id_type id) const;
        wt_t edge_weight(id_type from_id, id_type to_id) const;
        id_type insert_node();
        bool insert_edge(id_type from_id, id_type to_id, wt_t weight);
        bool delete_edge(id_type from_id, id_type to_id);
        wt_t set_edge_wt(id_type from_id, id_type to_id, wt_t weight);
        bool delete_node(id_type id);
        bool delete_all_edges();
        bool clear();
    };

    typedef unordered_map< id_type,pair<id_type,id_type> > id_node_map;
    typedef id_node_map::const_iterator ds_iterator;
    class disjoint_set {
    private:
        id_type st_num_sets;
        id_node_map inm_elems;
    public:
        disjoint_set();
        id_type size() const;
        id_type num_sets() const;
        ds_iterator begin() const;
        ds_iterator end() const;
        pair<ds_iterator,bool> make_set(id_type x);
        pair<ds_iterator,bool> find(id_type x);
        pair<ds_iterator,bool> join(id_type x, id_type y);
    };
    
    class binary_heap {
    private:
        vector< pair<id_type,wt_type> >vpiw_heap;
        unordered_map<id_type,id_type> umis_pos;
        bool b_max;
        bool compare(const pair<id_type,wt_type>& left,const pair<id_type,wt_type>& right) const;
    public:
        
        id_type size() const;
        bool empty() const;
        pair<id_type,wt_type> top() const;
        void update_key(const pair<id_type,wt_type>& piw_in);
        void insert(const pair<id_type,wt_type>& piw_in);
        void heapify_up(id_type pos);
        void heapify_down(id_type root);
        void pop();
        binary_heap(bool min);
        binary_heap(const vector<wt_type>& v,bool max);
    };
    
    
};

#endif	/* DATASTRUCTURES_H */

