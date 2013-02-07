/* 
 * File:   graph_properties.h
 * Author: bharath
 *
 * Created on May 15, 2012, 8:26 AM
 */

#ifndef GRAPH_PROPERTIES_H
#define	GRAPH_PROPERTIES_H
#include "typedefs.h"
#include "graph.h"
#include "statistics.h"
#include "random_graph.h"
using namespace std;


namespace CDLib 
{
    template <typename T>
    void get_degree_sequence(const graph& g,vector<T>& degrees, bool in_degrees);
    void get_degree_histogram(const graph& g,vector<id_type>& dist, bool in_degrees);
    double get_degree_distribution(const graph& g,vector<double>& dist,bool in_degrees);
    double get_degree_variance(const graph& g,bool in_degrees);
    void get_degree_assortativity_coefficient(const graph& g,bool in_degrees,vector<double>& assortativity);
    double get_degree_assortativity_coefficient(const graph& g,bool in_degrees);
    double get_rich_club_coefficient(const graph& g,id_type start_hub_degree);
    double normalized_rich_club_coefficient(const graph& g,id_type start_hub_degree_def);
    double get_poor_club_coefficient(const graph& g,id_type start_hub_degree);
    
    double kl_divergence_from_random_graph(const graph& g);
    double distance_from_random_graph(const graph& g, bool hellinger);  // Hellinger and Bhattacharyya Distance of Degree Distribution
    double connectivity_entropy(const graph& g);
    double graph_modularity(graph& g);
};

#endif	/* GRAPH_PROPERTIES_H */

