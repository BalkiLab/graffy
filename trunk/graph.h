/* 
 * File:   graph.h
 * Author: bharath
 *
 * Created on April 1, 2012, 12:00 PM
 */

#ifndef GRAPH_H
#define	GRAPH_H

#include "typedefs.h"
#include "datastructures.h"
#include "utility.h"

namespace CDLib {

    class graph {
    private:

        bool b_directed;
        bool b_weighted;
        bidirectional_label_map blm_labels;
        double_adjacency_map dam_backend;
        string graph_name;

    public:

        graph();
        graph(bool directed, bool weighted);

        bool is_directed() const;
        bool is_weighted() const;
        void set_graph_name(string gname);
        string get_graph_name() const;

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
        wt_t add_edge(const string& from_label, const string& to_label, wt_t weight);

        bool remove_node(id_type id);
        bool remove_node(const string& label);
        id_type remove_nodes(const node_set& nodes);
        id_type remove_nodes(const set<string>& nodes);
        id_type remove_isolates();

        bool remove_edge(id_type from_id, id_type to_id);
        wt_t remove_edge(const string& from_label, const string& to_label);
        id_type remove_adjacent_edges(id_type id);
        id_type remove_adjacent_edges(const string& label);
        id_type remove_edges(vector<pair<id_type, id_type> >& edges);
        id_type remove_edges(vector<pair<string, string> >& edges);

        bool set_edge_weight(id_type from_id, id_type to_id, wt_t weight);
        wt_t set_edge_weight(const string& from_label, const string& to_label, wt_t weight);

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

    struct edge {
        id_type from;
        id_type to;
        double weight;

        edge(id_type f, id_type t, double w) : from(f), to(t), weight(w) {
        }
    };

    class edge_iterator {
    private:
        const graph* g;
        id_type current_id;
        adjacent_edges_iterator current_edge;
    public:

        edge_iterator(const graph* gp) : g(gp) {
            current_id = 0;
            current_edge = g->out_edges_begin(current_id);
        }

        edge_iterator(const graph* gp, id_type cid, adjacent_edges_iterator cei) : g(gp), current_id(cid), current_edge(cei) {
        }

        edge_iterator& operator ++() {
            if (current_edge != g->out_edges_end(current_id)) {

                current_edge++;
                if (current_edge == g->out_edges_end(current_id)) {
                    if (current_id < g->get_num_nodes() - 1) {
                        current_id++;
                        current_edge = g->out_edges_begin(current_id);
                    } else current_edge = g->out_edges_end(current_id);
                }
            }
            return (*this);
        }

        edge_iterator& operator ++(int) {
            if (current_edge != g->out_edges_end(current_id)) {

                current_edge++;
                if (current_edge == g->out_edges_end(current_id)) {
                    if (current_id < g->get_num_nodes() - 1) {
                        current_id++;
                        current_edge = g->out_edges_begin(current_id);
                    } else current_edge = g->out_edges_end(current_id);
                }
            }
            return (*this);
        }

        edge operator *() {
            return edge(current_id, current_edge->first, current_edge->second);
        }

        bool operator !=(edge_iterator rhs) {
            return current_id != rhs.current_id && current_edge != rhs.current_edge;
        }
    };

    inline edge_iterator graph_edges_begin(const graph& g) {
        return edge_iterator(&g, 0, g.out_edges_begin(0));
    }

    inline edge_iterator graph_edges_end(const graph& g) {
        return edge_iterator(&g, g.get_num_nodes() - 1, g.out_edges_end(g.get_num_nodes() - 1));
    }

    class CSR {
    private:
        vector<id_type> xadj;
        vector<id_type> adjncy;
        vector<double> weights;
        double total_weight;
    public:
        typedef vector<id_type>::const_iterator edgeIter;

        inline id_type get_num_nodes() const {
            return xadj.size();
        }

        inline id_type get_num_edges() const {
            return adjncy.size();
        }

        inline double get_total_weight() const {
            if (weights.size()) return total_weight;
            else return 0;
        }

        inline id_type get_node_degree(id_type i) const {
            if (i < xadj.size() - 1) return xadj[i + 1] - xadj[i];
            else return (adjncy.size() - xadj[i]);
        }

        inline bool is_weighted() const {
            return weights.size();
        }

        inline edgeIter out_edges_begin(id_type id) const {
            if (id < xadj.size())return adjncy.begin() + xadj[id];
            else return adjncy.end();
        }

        inline edgeIter out_edges_end(id_type id) const {
            if (id + 1 < xadj.size())return adjncy.begin() + xadj[id + 1];
            else return adjncy.end();
        }

        inline edgeIter edges_begin() const {
            return adjncy.begin();
        }

        inline edgeIter edges_end() const {
            return adjncy.end();
        }

        edgeIter find_edge(id_type i, id_type j) const {
            if (i < xadj.size() && j < xadj.size()) {
                edgeIter it;
                if (i < xadj.size() - 1)it = lower_bound(adjncy.begin() + xadj[i], adjncy.begin() + xadj[i + 1], j);
                else it = lower_bound(adjncy.begin() + xadj[i], adjncy.end(), j);
                return it;
            }
            return adjncy.end();
        }

        inline double get_edge_weight(edgeIter eit) const {
            if (eit != adjncy.end() && weights.size())return *(weights.begin() + (eit - adjncy.begin()));
            else if (eit != adjncy.end()) return 1;
            return 0;
        }

        inline double get_edge_weight(id_type i, id_type j) const {
            return get_edge_weight(find_edge(i, j));
        }

        void populate_from_graph(const graph & g) {
            adjncy.clear();
            weights.clear();
            xadj.assign(g.get_num_nodes(), 0);
            id_type num_edges;
            total_weight = g.get_total_weight();
            if (g.is_directed()) num_edges = g.get_num_edges();
            else num_edges = (2 * g.get_num_edges()) - g.get_num_self_edges();
            adjncy.reserve(num_edges);
            if (g.is_weighted()) weights.reserve(num_edges);
            for (id_type i = 0; i < g.get_num_nodes(); i++) {
                if (i)xadj[i] = xadj[i - 1] + g.get_node_out_degree(i - 1);
                vector<id_type> out_edges;
                vector<double> out_edge_wts;
                out_edges.reserve(g.get_node_out_degree(i));
                for (adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i); aeit++)out_edges.push_back(aeit->first);
                sort(out_edges.begin(), out_edges.end());
                copy(out_edges.begin(), out_edges.end(), back_inserter(adjncy));
                if (g.is_weighted()) {
                    out_edge_wts.reserve(g.get_node_out_degree(i));
                    for (id_type j = 0; j < out_edge_wts.size(); j++)out_edge_wts.push_back(g.get_edge_weight(i, j));
                    copy(out_edge_wts.begin(), out_edge_wts.end(), back_inserter(weights));
                }
            }
        }

        CSR() {
        }

        CSR(const graph & g) {
            populate_from_graph(g);
        }
    };

};

#endif	/* GRAPH_H */

