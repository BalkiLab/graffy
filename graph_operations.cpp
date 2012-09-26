/* 
 * File:   graph_operations.cpp
 * Author: bharath
 * 
 * Created on April 18, 2012, 6:03 AM
 */

#include "graph_operations.h"
using namespace CDLib;


id_type CDLib::extract_subgraph(const graph& g, node_set& nodes, graph& sg)
{
    sg.clear();
    if(sg.is_directed()!=g.is_directed() && sg.is_weighted()!=g.is_weighted()) return 0;
    for(node_set::iterator nit = nodes.begin(); nit != nodes.end(); nit++)
    {
        string first_label = g.get_node_label(*nit);
        sg.add_node(first_label);
        for(adjacent_edges_iterator aeit = g.out_edges_begin(*nit);aeit != g.out_edges_end(*nit); aeit++)
        {
            if(nodes.find(aeit->first)!=nodes.end())
            {
                string second_label =  g.get_node_label(aeit->first);
                sg.add_node(second_label);
                sg.add_edge(first_label,second_label,aeit->second);
            }
        }
    }
    return sg.get_num_edges();
}

id_type CDLib::copy_graph(const graph& src, graph& dst)
{
    // Naive implementation to Copy a graph. May be inefficient
    // Returns the number of nodes in the new graph
    dst.clear();
    if(src.is_directed()!=dst.is_directed() && src.is_weighted()!=dst.is_weighted()) return 0;
    for(id_type i=0;i < src.get_num_nodes();i++){
        dst.add_node(src.get_node_label(i));
    }
    for(id_type i=0;i < src.get_num_nodes();i++){
        for(adjacent_edges_iterator aeit = src.out_edges_begin(i);aeit != src.out_edges_end(i); aeit++){
            dst.add_edge(src.get_node_label(i),src.get_node_label(aeit->first),aeit->second);
        }
    }
    return dst.get_num_nodes();
}

void CDLib::sample_graph(const graph&g,node_set& seeds,id_type hop_dist,graph& sample)
{
    node_set new_nodes;
    for(node_set::iterator nit = seeds.begin(); nit != seeds.end();nit++)
    {
        vector<double> distances(g.get_num_nodes(),numeric_limits<double>::infinity());
        queue<id_type> bfs_queue;
        bfs_queue.push(*nit);
        distances[*nit] = 0;
        new_nodes.insert(*nit);
        while(!bfs_queue.empty())
        {
            id_type curr = bfs_queue.front();
            bfs_queue.pop();
            for(adjacent_edges_iterator aeit = g.out_edges_begin(curr);aeit != g.out_edges_end(curr);aeit++)
            {
                if(distances[aeit->first] == numeric_limits<double>::infinity() && distances[curr] <hop_dist) 
                {
                    distances[aeit->first] = distances[curr] + 1;
                    new_nodes.insert(aeit->first);
                    bfs_queue.push(aeit->first);
                }
            }
        }
    }
    extract_subgraph(g,new_nodes,sample);
}

void CDLib::multiply_vector_transform(const graph& g,vector<double>& invec,double (*wt_transform_func)(const graph&g,id_type,id_type,double),vector<double>& outvec)
{
    if(invec.size()==g.get_num_nodes())
    {
        outvec.assign(g.get_num_nodes(),0);
        for(id_type i=0;i<g.get_num_nodes();i++)
            for(adjacent_edges_iterator aeit = g.in_edges_begin(i);aeit != g.in_edges_end(i);aeit++)
                outvec[i] += invec[aeit->first]*wt_transform_func(g,i,aeit->first,aeit->second);
    }
}

void CDLib::random_walk(const graph& g,vector<double>& invec,id_type t,double (*wt_transform_func)(const graph&g,id_type,id_type,double),vector<double>& outvec)
{
    if(invec.size()==g.get_num_nodes() && t>=0)
    {
        vector<double> qtemp(invec);
        for(id_type i=0;i<t;i++)
        {
            multiply_vector_transform(g,qtemp,wt_transform_func,outvec);
            qtemp = outvec;
        }
    }
}

double CDLib::transform_func_nop(const graph& g, id_type i,id_type j, double wt) { return wt; }
double CDLib::transform_func_row_stochastic(const graph& g, id_type i,id_type j, double wt) { return wt/g.get_node_out_weight(i); }
double CDLib::transform_func_column_stochastic(const graph& g, id_type i,id_type j, double wt) { return wt/g.get_node_in_weight(j); }
double CDLib::transform_func_max_rowcol_stochastic(const graph& g, id_type i,id_type j, double wt) { return wt/max(g.get_node_out_weight(i),g.get_node_in_weight(j)); }
double CDLib::transform_func_min_rowcol(const graph& g, id_type i,id_type j, double wt) { return wt/min(g.get_node_out_weight(i),g.get_node_in_weight(j)); }
double CDLib::transform_func_laplacian(const graph& g, id_type i,id_type j, double wt) 
{ 
    if(i == j) return g.get_node_out_degree(i) - wt;
    else return -wt;
}
double CDLib::transform_func_normalized_laplacian(const graph& g, id_type i,id_type j, double wt) 
{ 
    return 1 - wt/g.get_node_out_degree(i);
}
double CDLib::transform_func_modularity(const graph& g, id_type i,id_type j, double wt)
{
    return wt - ((g.get_node_out_degree(i)*g.get_node_in_degree(j))/(2*g.get_num_edges()));
}