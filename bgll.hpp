/* 
 * File:   bgll.hpp
 * Author: bharath
 *
 * Created on 25 January, 2013, 5:19 PM
 */

#ifndef BGLL_HPP
#define	BGLL_HPP


#include "graph.h"
#include "community_tools.h"

namespace CDLib {

    class bgll_objective {
    public:
        virtual void init(const graph& orig_graph, const vector<vector<id_type> >& hiercomms, const graph& curr_graph, const vector<id_type>& curr_comms) = 0;
        virtual double objval(const graph& curr_graph, const vector<id_type>& curr_comms) const = 0;
        virtual double compute_gain(const graph& curr_graph, const vector<id_type>& labels, id_type vertex, id_type dst_comm) const = 0;
        virtual void detach_node(const graph& curr_graph, const vector<id_type>& labels, id_type vertex) = 0;
        virtual void attach_node(const graph& curr_graph, const vector<id_type>& labels, id_type vertex, id_type dst_comm) = 0;
    };

    id_type bgll_find_best_community_for_node(const graph& curr_graph, const vector<id_type>& labels, const bgll_objective& book, id_type vertex) {
        double max_node_gain = -numeric_limits<double>::infinity();
        id_type node_max_gain_comm_index = labels[vertex];
        for (adjacent_edges_iterator aeit = curr_graph.out_edges_begin(vertex); aeit != curr_graph.out_edges_end(vertex); aeit++) {
            if (labels[vertex] != labels[aeit->first]) {
                double curr_node_gain = book.compute_gain(curr_graph, labels, vertex, labels[aeit->first]);
                if (curr_node_gain > max_node_gain) {
                    max_node_gain = curr_node_gain;
                    node_max_gain_comm_index = labels[aeit->first];
                }
            }
        }
        //cout << node_max_gain_comm_index << " " << max_node_gain <<endl;
        return node_max_gain_comm_index;
    }

    id_type bgll_vertex_mover_single_pass(const graph& curr_graph, vector<id_type>& labels, bgll_objective& book) {
        id_type num_nodes_moved = 0;
        for (id_type i = 0; i < curr_graph.get_num_nodes(); i++) {
            book.detach_node(curr_graph, labels, i); //Detach node from its community
            //cout << labels[i] << " ";
            id_type new_node_community = bgll_find_best_community_for_node(curr_graph, labels, book, i);
            if (new_node_community != labels[i]) {
                num_nodes_moved++;
                //Insert node into new community
                book.attach_node(curr_graph, labels, i, new_node_community);
                labels[i] = new_node_community;
            }
        }
        return num_nodes_moved;
    }

    double bgll_vertex_mover_optimizer(const graph& g, vector<id_type>& labels, bgll_objective& book) {
        double start_obj_val = book.objval(g, labels), curr_objval = -numeric_limits<double>::infinity();
        id_type num_nodes_moved = 0;
        do {
            vector<id_type> temp_comms(labels);
            curr_objval = book.objval(g, temp_comms);
            num_nodes_moved += bgll_vertex_mover_single_pass(g, temp_comms, book);
            //Rollback communities and exit if modularity decreases
            if (curr_objval < start_obj_val) return start_obj_val;
            //Commit changes, book should have already been updated
            copy(temp_comms.begin(), temp_comms.end(), labels.begin());
        } while (num_nodes_moved && curr_objval < start_obj_val);
        return curr_objval;
    }

    void bgll_recover_communities(vector<id_type>& labels, vector< vector<id_type> >& hier_comms) {
        vector<node_set> last_comms;
        convert_labels_to_communities(hier_comms.back(), last_comms);
        reindex_communities(labels);
        hier_comms.push_back(vector<id_type > (hier_comms.back().size(), 0));
        for(id_type i=0;i<labels.size();i++){
            for(node_set::iterator nit=last_comms[i].begin();nit!=last_comms[i].end();nit++){
                hier_comms[hier_comms.size()-2][*nit] = labels[i];
            }
        }
    }

    void bgll_collapse_nodes(graph& orig_graph, vector<id_type>& labels) {
        graph g(0,1);
        for (id_type i = 0; i < labels.size(); i++)
            while (g.get_num_nodes() <= labels[i]) g.add_node();
        for (id_type i = 0; i < orig_graph.get_num_nodes(); i++) {
            id_type node_comm = labels[i];
            for (adjacent_edges_iterator aeit = orig_graph.out_edges_begin(i); aeit != orig_graph.out_edges_end(i); aeit++) {
                id_type neigh_comm = labels[aeit->first];
                double new_edge_weight = g.get_edge_weight(node_comm, neigh_comm) + (aeit->second/2);
                g.set_edge_weight(node_comm, neigh_comm,new_edge_weight);
            }
        }
        cout <<orig_graph.get_num_nodes() << " " << orig_graph.get_total_weight() << " " << g.get_num_nodes() << " "<< g.get_total_weight() <<endl;
        orig_graph.clear();
        copy_graph(g,orig_graph);
        labels.assign(orig_graph.get_num_nodes(), 0);
        for (id_type i = 0; i < orig_graph.get_num_nodes(); i++)
            labels[i] = i;
    }

    void bgll_community_detection(const graph&g, vector<id_type>& init_comms, vector< vector<id_type> >& hier_comms, bgll_objective& book) {
        hier_comms.clear();
        //Initialize the initial_partition
        vector<id_type> curr_comms(g.get_num_nodes(), 0);
        if (init_comms.size() == g.get_num_nodes()) copy(init_comms.begin(), init_comms.end(), curr_comms.begin());
        else for (id_type i = 0; i < g.get_num_nodes(); i++) curr_comms[i] = i;
        hier_comms.push_back(curr_comms);
        //Create a copy of the graph
        graph curr_graph(g);
        curr_graph.convert_to_weighted();
        double orig_obj_val = book.objval(curr_graph, curr_comms);
        //The Main Loop
        while (1) {
            id_type curr_graph_size = curr_graph.get_num_nodes();
            //Initialize the book for the new level
            book.init(g, hier_comms, curr_graph, curr_comms);
            //Run the vertex mover optimization
            double new_obj_val = bgll_vertex_mover_optimizer(curr_graph, curr_comms, book);
            //Reindex the community memberships
            bgll_recover_communities(curr_comms, hier_comms);
            //Shrink the graph
            bgll_collapse_nodes(curr_graph, curr_comms);
            if (curr_graph.get_num_nodes() == curr_graph_size || new_obj_val < orig_obj_val) break;
        }
    }

    class bgll_modularity : public bgll_objective {
    private:
        unordered_map<id_type, double> internal_weight;
        unordered_map<id_type, double> total_weight;
        unordered_map<id_type, double> node_link_weights;
        double resolution_param;
        double self_weight;
    public:
        bgll_modularity () :resolution_param(1) {} 
        bgll_modularity(double resolution_p) : resolution_param(resolution_p){ }
        virtual void init(const graph& orig_graph, const vector<vector<id_type> >& hiercomms, const graph& curr_graph, const vector<id_type>& curr_comms) {
            for (id_type i = 0; i < curr_graph.get_num_nodes(); i++) {
                id_type node_comm = curr_comms[i];
                for (adjacent_edges_iterator aeit = curr_graph.out_edges_begin(i); aeit != curr_graph.out_edges_end(i); aeit++) {
                    id_type neigh_comm = curr_comms[aeit->first];
                    if (node_comm == neigh_comm) {
                        pair < unordered_map<id_type, double>::iterator, bool> ret = internal_weight.insert(make_pair(node_comm, 0));
                        ret.first->second += aeit->second;
                    }
                }
                pair < unordered_map<id_type, double>::iterator, bool> ret = internal_weight.insert(make_pair(node_comm, 0));
                ret.first->second += curr_graph.get_node_out_weight(i);
            }
        }

        virtual double objval(const graph& curr_graph, const vector<id_type>& curr_comms) const {
            double modularity = 0;
            for (unordered_map<id_type, double>::const_iterator it = total_weight.begin(); it != total_weight.end(); it++) {
                double int_weight = 0;
                unordered_map<id_type, double>::const_iterator cit = internal_weight.find(it->first);
                if (cit != internal_weight.end()) int_weight = cit->second;
                modularity += ((resolution_param * int_weight) - pow(it->second / (2 * curr_graph.get_total_weight()), 2));
            }
            return modularity;
        }

        virtual double compute_gain(const graph& curr_graph, const vector<id_type>& labels, id_type vertex, id_type dst_comm) const {
            unordered_map<id_type, double>::const_iterator it = node_link_weights.find(dst_comm);
            double gain = 0;
            if (it != node_link_weights.end())gain += resolution_param * 2 * it->second;
            it = total_weight.find(dst_comm);
            if (it != total_weight.end()) gain -= ((curr_graph.get_node_out_weight(vertex) * it->second) / curr_graph.get_total_weight());
            return gain;
        }

        virtual void detach_node(const graph& curr_graph, const vector<id_type>& labels, id_type vertex) {
            node_link_weights.clear();
            self_weight = 0;
            double weight_inside = 0;
            for (adjacent_edges_iterator aeit = curr_graph.out_edges_begin(vertex); aeit != curr_graph.out_edges_end(vertex); aeit++) {
                if (labels[aeit->first] != labels[vertex]) {
                    pair < unordered_map<id_type, double>::iterator, bool> ret = node_link_weights.insert(make_pair(labels[aeit->first], 0));
                    ret.first->second += aeit->second;
                } else if (vertex != aeit->first) {
                    weight_inside += aeit->second;
                } else self_weight = aeit->second;
            }
            unordered_map<id_type, double>::iterator it1 = total_weight.find(labels[vertex]), it2 = internal_weight.find(labels[vertex]);
            if (it1 != total_weight.end() && it2 != total_weight.end()) {
                it1->second -= curr_graph.get_node_out_weight(vertex);
                it2->second -= (2 * weight_inside + self_weight);
            }
        }

        virtual void attach_node(const graph& curr_graph, const vector<id_type>& labels, id_type vertex, id_type dst_comm) {
            if (labels[vertex] != dst_comm) {
                unordered_map<id_type, double>::iterator it1 = total_weight.find(dst_comm), it2 = internal_weight.find(dst_comm), it3 = node_link_weights.find(dst_comm);
                if (it1 != total_weight.end() && it2 != total_weight.end() && it3 != node_link_weights.end()) {
                    it1->second += curr_graph.get_node_out_weight(vertex);
                    it2->second += (2 * it3->second + self_weight);
                }
            }
        }
    };

};

#endif	/* BGLL_HPP */

