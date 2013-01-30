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

double CDLib::get_degree_distribution(const graph& g,vector<double>& dist, bool in_degrees)
{
    vector<id_type> degree_hist;
    dist.clear();
    double expectation = 0;
    get_degree_histogram(g,degree_hist,in_degrees);
    dist.assign(degree_hist.size(),0);
    for (id_type i=0; i<degree_hist.size(); i++) {
        dist[i] = (double)degree_hist[i]/(double)g.get_num_nodes();
        expectation += i * dist[i];
    }
    return expectation;
}

double CDLib::get_degree_variance(const graph& g,bool in_degrees)
{
    if (g.get_num_nodes() <= 0)
        return 0;
    vector<id_type> degrees;
    get_degree_sequence(g,degrees,in_degrees);
    return variance(degrees);
}

void CDLib::get_degree_assortativity_coefficient(const graph& g,bool in_degrees,vector<double>& assortativity)
{
/* Gives the assortativity contributions from each node */    
/* The calculation is based on Classifying Complex Networks using Unbiased Local Assortativity, 2010 paper by Piraveenan et al. */        
    assortativity.clear();
    vector<double> degree_distt;
    double mean = get_degree_distribution(g,degree_distt,in_degrees);
    vector<double> excess_degree_distt(degree_distt.size()-1,0);
    double edd_mean = 0,edd_variance = 0;
    for (id_type i=0;i<excess_degree_distt.size();i++) {
        excess_degree_distt[i] = ((i+1) * degree_distt[i+1])/mean;
        double imean = i*excess_degree_distt[i];
        edd_mean += imean;
        edd_variance += i*imean;
    }
    edd_variance -= (edd_mean*edd_mean);
    if (edd_variance == 0) {        
        assortativity.assign(g.get_num_nodes(),(1/g.get_num_nodes()));
        return;
    }
    else {      
        assortativity.assign(g.get_num_nodes(),0);    
    }
    id_type edge_factor = 2 * g.get_num_edges();
#pragma omp parallel for schedule(dynamic,20) shared(g,assortativity)        
    for(id_type i=0;i<g.get_num_nodes();i++) {
        double avg_excess_degree_neighbour = 0;
        if (in_degrees) {
            for(adjacent_edges_iterator aeit = g.in_edges_begin(i);aeit != g.in_edges_end(i);aeit++) {
                avg_excess_degree_neighbour += g.get_node_in_degree(aeit->first);
            }
            avg_excess_degree_neighbour = (avg_excess_degree_neighbour/g.get_node_in_degree(i))-1;
            assortativity[i] = ((g.get_node_in_degree(i) - 1) * g.get_node_in_degree(i) * (avg_excess_degree_neighbour - edd_mean))/(edd_variance * edge_factor);
        }
        else {
            for(adjacent_edges_iterator aeit = g.out_edges_begin(i);aeit != g.out_edges_end(i);aeit++) {
                avg_excess_degree_neighbour += g.get_node_out_degree(aeit->first);
            }
            avg_excess_degree_neighbour = (avg_excess_degree_neighbour/g.get_node_out_degree(i))-1;
            assortativity[i] = ((g.get_node_out_degree(i) - 1) * g.get_node_out_degree(i) * (avg_excess_degree_neighbour - edd_mean))/(edd_variance * edge_factor);
        }
    }
}

double CDLib::get_degree_assortativity_coefficient(const graph& g,bool in_degrees)
{
/* The calculation is based on Assortative Mixing in Networks, 2002 paper by MEJ Newman */        
    if (g.get_num_nodes() <= 0)
        return 0;
    id_type edge_factor = 2 * g.get_num_edges();
    double sum=0,prod_sum=0,square_sum=0;
#pragma omp parallel for schedule(dynamic,20) shared(g) reduction(+:sum,prod_sum,square_sum)        
    for(id_type i=0;i<g.get_num_nodes();i++) {
        if (in_degrees) {
            for(adjacent_edges_iterator aeit = g.in_edges_begin(i);aeit != g.in_edges_end(i);aeit++) {
                sum += g.get_node_in_degree(i) + g.get_node_in_degree(aeit->first);
                prod_sum += g.get_node_in_degree(i) * g.get_node_in_degree(aeit->first);
                square_sum += pow(g.get_node_in_degree(i),2) + pow(g.get_node_in_degree(aeit->first),2);
            }
        }
        else {
            for(adjacent_edges_iterator aeit = g.out_edges_begin(i);aeit != g.out_edges_end(i);aeit++) {
                sum += g.get_node_out_degree(i) + g.get_node_out_degree(aeit->first);
                prod_sum += g.get_node_out_degree(i) * g.get_node_out_degree(aeit->first);
                square_sum += pow(g.get_node_out_degree(i),2) + pow(g.get_node_out_degree(aeit->first),2);
            }
        }
    }
    sum /= edge_factor;             square_sum /= edge_factor;      prod_sum /= g.get_num_edges();
    double tmp_sq = sum * sum;
    return (prod_sum - tmp_sq)/(square_sum - tmp_sq);
}

double CDLib::get_rich_club_coefficient(const graph& g,id_type hub_degree_def)
{
/* This implementation is accordance to Detecting rich-club ordering in complex networks, 2006 paper by Colizza et al. */    
    id_type num_rich_nodes = 0, num_rich_edges = 0;
#pragma omp parallel for schedule(dynamic,20) shared(g) reduction(+:num_rich_nodes,num_rich_edges)         
    for(id_type i=0;i<g.get_num_nodes();i++) {
        if (g.get_node_out_degree(i) > hub_degree_def) {
            num_rich_nodes++;
            for(adjacent_edges_iterator aeit = g.out_edges_begin(i);aeit != g.out_edges_end(i);aeit++) {
                if (g.get_node_out_degree(aeit->first) > hub_degree_def) num_rich_edges++;
            }
        }   
    }
    return (double)(2 * num_rich_edges)/(num_rich_nodes * (num_rich_nodes-1));
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
    vector<double> distt1, distt2;
    double lambda = (2 * g.get_num_edges())/g.get_num_nodes();
    id_type degree = 0;
    double probability = exp(-1*lambda);
    double factorial = 1;
    while (((probability != HUGE_VAL) || (probability != -HUGE_VAL)) && ((degree < lambda) || (probability > 0.001))) {
        distt2.push_back(probability);
        degree++;
        factorial *= degree;
        probability = (pow(lambda,degree) * exp(-1*lambda))/factorial;
    }
    get_degree_distribution(g,distt1,0);
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
