/* 
* File:   double_adjacency_map.h
* Author: bharath
*
* Created on April 1, 2012, 1:50 AM
*/

#ifndef DOUBLE_ADJACENCY_MAP_H
#define	DOUBLE_ADJACENCY_MAP_H
#include "typedefs.h"
using namespace std;
namespace CDLib
{
       typedef unordered_map<id_type,wt_t> adjacent_edge_sequence;
       typedef vector<adjacent_edge_sequence> adjacency_map;
       typedef typename adjacent_edge_sequence::const_iterator adjacent_edges_iterator;
       class double_adjacency_map 
       {
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
           bool insert_edge(id_type from_id, id_type to_id,wt_t weight);
           bool delete_edge(id_type from_id,id_type to_id);
           wt_t set_edge_wt(id_type from_id,id_type to_id, wt_t weight);
           bool delete_node(id_type id);
           bool delete_all_edges();
           bool clear();
       };
};

#endif	/* DOUBLE_ADJACENCY_MAP_H */

