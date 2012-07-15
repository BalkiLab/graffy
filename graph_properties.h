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
using namespace std;


namespace CDLib 
{
    string get_graph_details(const graph& g);
    void get_degree_distribution(const graph& g,vector<id_type>& dist,bool in_degrees);
};

#endif	/* GRAPH_PROPERTIES_H */

