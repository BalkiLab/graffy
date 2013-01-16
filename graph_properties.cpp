/* 
 * File:   graph_properties.cpp
 * Author: bharath
 * 
 * Created on May 15, 2012, 8:26 AM
 */

#include "graph_properties.h"
#include "divisive_algorithms.h"

void CDLib::get_degree_sequence(const graph& g,vector<id_type>& degrees, bool in_degrees)
{
    degrees.clear();
    degrees.assign(g.get_num_nodes(),0);
    for(id_type i=0;i<g.get_num_nodes();i++)
        degrees[i] = (in_degrees)? g.get_node_in_degree(i) : g.get_node_out_degree(i);
}

void CDLib::get_degree_histogram(const graph& g,vector<id_type>& dist, bool in_degrees)
{
    vector<id_type> degrees;
    get_degree_sequence(g,degrees,in_degrees);
    get_discrete_distribution<id_type>(degrees,dist);
}

void CDLib::get_degree_distribution(const graph& g,vector<double>& dist, bool in_degrees)
{
    vector<id_type> degree_hist;
    dist.clear();
    get_degree_histogram(g,degree_hist,in_degrees);
    dist.assign(degree_hist.size(),0);
    for (id_type i=0; i<degree_hist.size(); i++)
        dist[i] = (double)degree_hist[i]/(double)g.get_num_nodes();
}

double CDLib::get_degree_variance(const graph& g,bool in_degrees)
{
    if (g.get_num_nodes() <= 0)
        return 0;
    vector<id_type> degrees;
    get_degree_sequence(g,degrees,in_degrees);
    return variance(degrees);
}

double CDLib::get_degree_assortativity_coefficient(const graph& g,bool in_degrees)
{
    if (g.get_num_nodes() <= 0)
        return 0;
    vector<id_type> degrees;
    get_degree_sequence(g,degrees,in_degrees);
    double assortative_cov = 0,max_assortative_cov = 0;
    id_type edge_factor = 2 * g.get_num_edges();
#pragma omp parallel for schedule(dynamic,10) shared(g,degrees,edge_factor) reduction(+:assortative_cov,max_assortative_cov)     
    for(id_type i=0;i<g.get_num_nodes();i++) {
        for(id_type j=0;j<g.get_num_nodes();j++) {
            double deg_prod = degrees[i] * degrees[j];
            assortative_cov += (g.get_edge_weight(i,j) - (deg_prod/edge_factor)) * deg_prod;
            max_assortative_cov += ((degrees[i] * kronecker_delta(degrees[i],degrees[j])) - (deg_prod/edge_factor)) * deg_prod;
        }
    }
    return (assortative_cov/max_assortative_cov);
}

double CDLib::kl_divergence_from_random_graph(const graph& g)
{
//    Reports symmetric KL-Divergence of the Degree Distribution from Random graph of same size. 
    graph er_graph(0,0);
    generate_erdos_renyi_graph(er_graph, g.get_num_nodes(),g.get_num_edges());
    vector<double> distt1, distt2;
    get_degree_distribution(g,distt1,0);
    get_degree_distribution(er_graph,distt2,0);
    return kl_divergence_symmetric(distt1,distt2);
}
double CDLib::distance_from_random_graph(const graph& g, bool hellinger)
{
//    If bool is set, it return Hellinger Distance of the Degree Distribution from Random Graph of same size,
//    else it returns Bhattacharyya Distance of the Degree Distribution.
    graph er_graph(0,0);
    generate_erdos_renyi_graph(er_graph, g.get_num_nodes(),g.get_num_edges());
    vector<double> distt1, distt2;
    get_degree_distribution(g,distt1,0);
    get_degree_distribution(er_graph,distt2,0);
    if (hellinger)
        return hellinger_distance(distt1,distt2);
    else
        return bhattacharyya_distance(distt1,distt2);
}

double CDLib::connectivity_entropy(const graph& g)
{
    // Connective Entropy of the Network or Information Entropy of the Network
    double entropy = 0;
    if ((2 * g.get_num_edges()) > 0){
        for (id_type i = 0; i < g.get_num_nodes(); i++){
            double prob = (double)g.get_node_out_degree(i)/(2 * g.get_num_edges());
            if (prob != 0){ entropy += prob * (log(prob)/log(2)); }
        }
        entropy *= -1;
    }
    entropy /= log(g.get_num_nodes())/log(2);   // Normalization of Entropy
    return entropy;
}

double CDLib::graph_modularity(graph& g)
{
    /* Functions returns the modularity of the whole graph considering the whole graph as a single community */
    double denominator = 2 * g.get_num_edges();
    double sum = 0;
    if (g.is_directed()){
        for (id_type i=0; i<g.get_num_nodes();i++){
            for (id_type j=0; j<g.get_num_nodes();j++){
                sum += g.get_edge_weight(i,j) - ((g.get_node_out_degree(i) * g.get_node_out_degree(j))/denominator);
            }
        } 
    }
    else{
        for (id_type i=0; i<g.get_num_nodes();i++){
            for (id_type j=i+1; j<g.get_num_nodes();j++){
                sum += 2 * (g.get_edge_weight(i,j) - ((g.get_node_out_degree(i) * g.get_node_out_degree(j))/denominator));
            }
        }
        for (id_type i=0; i<g.get_num_nodes();i++){
            sum += g.get_edge_weight(i,i) - ((g.get_node_out_degree(i) * g.get_node_out_degree(i))/denominator);
        }
    }
    return sum;
}
