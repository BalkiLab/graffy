/* 
 * File:   CSR.hpp
 * Author: bharath
 *
 * Created on August 17, 2012, 11:00 AM
 */

#ifndef CSR_HPP
#define	CSR_HPP

#include "graph.h"

class CSR {
private:
    vector<id_type> xadj;
    vector<id_type> adjncy;
    vector<double> weights;
public:
    typedef vector<id_type>::const_iterator edgeIter;

    inline id_type get_num_nodes() const {
        return xadj.size();
    }

    inline id_type get_num_edges() const {
        return adjncy.size();
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

    inline double get_edge_weight(edgeIter & eit) const {
        if (eit != adjncy.end() && weights.size())return weights.begin() + eit - adjncy.begin();
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
        if (g.is_directed()) num_edges = g.get_num_edges();
        else num_edges = (2 * g.get_num_edges()) - g.get_num_self_edges();
        adjncy.reserve(num_edges);
        if (g.is_weighted()) weights.reserve(num_edges);
        for (id_type i = 0; i < g.get_num_nodes(); i++) {
            if (i)xadj[i] = xadj[i - 1] + g.get_node_out_degree(i - 1);
            vector<id_type> out_edges;
            vector<double> out_edge_wts;
            out_edges.reserve(g.get_node_out_degree(i));
            for (adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(); aeit++)out_edges.push_back(aeit->first);
            sort(out_edges.begin(), out_edges.end());
            copy(out_edges.begin(), out_edges.end(), back_inserter(adjncy));
            if (g.is_weighted()) {
                out_edge_wts.reserve(g.get_node_out_degree(i));
                for (id_type j = 0; j < out_edge_wts.size(); j++)out_edge_wts.push_back(g.get_edge_weight(i, j));
                copy(out_edge_wts.begin(), out_edge_wts.end(), back_inserter(weights));
            }
        }
    }

    CSR(const graph & g) {
        populate_from_graph(g);
    }
};


#endif	/* CSR_HPP */

