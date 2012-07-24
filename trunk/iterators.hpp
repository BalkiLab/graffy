/* 
 * File:   iterators.hpp
 * Author: bharath
 *
 * Created on 24 July, 2012, 12:35 PM
 */

#ifndef ITERATORS_HPP
#define	ITERATORS_HPP

#include "graph.h"
using namespace std;

namespace CDLib
{
    struct edge
    {
        id_type from;
        id_type to;
        double weight;
        edge(id_type f, id_type t, double w) : from(f),to(t),weight(w) {}
    };
    class edge_iterator
    {
    private:
        const graph* g;
        id_type current_id;
        adjacent_edges_iterator current_edge;
    public:
        edge_iterator(const graph* gp) : g(gp)
        {
            current_id = 0;
            current_edge = g->out_edges_begin(current_id);
        }
        edge_iterator(const graph* gp,id_type cid,adjacent_edges_iterator cei) : g(gp),current_id(cid),current_edge(cei) {}
        edge_iterator& operator ++()
        {
            if(current_edge != g->out_edges_end(current_id))
            {
                
                current_edge++;
                if(current_edge == g->out_edges_end(current_id))
                {
                    if(current_id < g->get_num_nodes()-1)
                    {
                        current_id++;
                        current_edge = g->out_edges_begin(current_id);
                    }
                    else current_edge = g->out_edges_end(current_id);
                }
            }
            return (*this);
        }
        edge_iterator& operator ++(int)
        {
            if(current_edge != g->out_edges_end(current_id))
            {
                
                current_edge++;
                if(current_edge == g->out_edges_end(current_id))
                {
                    if(current_id < g->get_num_nodes()-1)
                    {
                        current_id++;
                        current_edge = g->out_edges_begin(current_id);
                    }
                    else current_edge = g->out_edges_end(current_id);
                }
            }
            return (*this);
        }
        edge operator *()
        {
            return edge(current_id,current_edge->first,current_edge->second);
        }
    };
}

edge_iterator graph_edges_begin(const graph& g)
{
    return edge_iterator(&g,0,g.out_edges_begin(0));
}

edge_iterator graph_edges_end(const graph& g)
{
    return edge_iterator(&g,g.get_num_nodes()-1,g.out_edges_end(g.get_num_nodes()-1));
}
#endif	/* ITERATORS_HPP */

