/* 
 * File:   graph_properties.cpp
 * Author: bharath
 * 
 * Created on May 15, 2012, 8:26 AM
 */

#include "graph_properties.h"


string CDLib::get_graph_details(const graph& g)
{
    ostringstream oss;
    oss << "Number of Nodes: " << g.get_num_nodes() << endl;
    oss << "Number of Edges: " << g.get_num_edges() << endl;
    oss << "Number of Self Edges: " << g.get_num_self_edges() << endl;
    oss << "Total Weight: " << g.get_total_weight() << endl;
    oss << "Weight of Self Edges: " << g.get_self_edges_weight() << endl;
    vector<id_type> in_degrees(g.get_num_nodes(),0),out_degrees(g.get_num_nodes(),0);
    vector<double> in_weights(g.get_num_nodes(),0.0),out_weights(g.get_num_nodes(),0.0);
    for(id_type i=0;i<g.get_num_nodes();i++)
    {
        in_degrees[i] = g.get_node_in_degree(i);
        out_degrees[i] = g.get_node_out_degree(i);
        in_weights[i] = g.get_node_in_weight(i);
        out_weights[i] = g.get_node_out_weight(i);
    }
    oss << "In Degrees" << endl;
    oss << "Min\tMedian\tMax\tMean\tVariance" << endl;
    oss << statistics_string(in_degrees,"\t") << endl;
    oss << "Out Degrees" << endl;
    oss << "Min\tMedian\tMax\tMean\tVariance" << endl;
    oss << statistics_string(out_degrees,"\t")  << endl;
    oss << "In Weights" << endl;
    oss << "Min\tMedian\tMax\tMean\tVariance" << endl;
    oss << statistics_string(in_weights,"\t") << endl;
    oss << "Out Weights" << endl;
    oss << "Min\tMedian\tMax\tMean\tVariance" << endl;
    oss << statistics_string(out_weights,"\t")  << endl;
    return oss.str();
}

void CDLib::get_degree_distribution(const graph& g,vector<id_type>& dist, bool in_degrees)
{
    vector<id_type> degrees(g.get_num_nodes(),0);
    for(id_type i=0;i<g.get_num_nodes();i++)
        degrees[i] = (in_degrees)? g.get_node_in_degree(i) : g.get_node_out_degree(i);
    get_discrete_distribution<id_type>(degrees,dist);
}
