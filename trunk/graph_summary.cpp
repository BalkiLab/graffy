/* 
 * File:   graph_summary.cpp
 * Author: sudip
 * 
 * Created on 18 October, 2012, 00:51 PM
 */

#include "graph_summary.h"


using namespace CDLib;

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

string CDLib::get_graph_details_components(const graph& g)
{
    ostringstream oss;
    vector<node_set> components;
    id_type num_comp = get_connected_components_undirected(g,components);
    oss << "Number of Components : " << num_comp << endl;
    oss << "Average Component Size : " << ((double)g.get_num_nodes()/num_comp) << endl;
    node_set members;
    id_type large_size = get_largest_connected_component(g,members);
    oss << "Size of Largest Connected Component : " << large_size << endl;
    oss << "Fraction of Nodes in LCC : " << ((double)large_size/g.get_num_nodes()) << endl;
    return oss.str();
}

string CDLib::get_graph_details_efficiency(const graph& g)
{
    ostringstream oss;
    oss << "Diameter : " << diameter(g) << endl;
    oss << "Small-World Efficiency : " << efficiency_sw_global(g,0) << endl;
    oss << "Connectivity Entropy (i.e. Entropy based on Degree Distt.) : " << connectivity_entropy(g) << endl;
    oss << "Path Entropy (i.e. Entropy based on Distance Distt.) : " << path_entropy(g) << endl;
    return oss.str();
}

string CDLib::get_graph_comparisons(const graph& g)
{
    ostringstream oss;
    oss << "KL-Divergence of Degree Distribution from Random Graph : " << kl_divergence_from_random_graph(g) << endl;
    oss << "Hellinger Distance of Degree Distribution from Random Graph : " << distance_from_random_graph(g,1) << endl;
    return oss.str();
}
