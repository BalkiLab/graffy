/* 
 * File:   graph_summary.h
 * Author: sudip
 *
 * Created on 18 October, 2012, 12:56 AM
 */

#ifndef GRAPH_SUMMARY_H
#define	GRAPH_SUMMARY_H

#include "typedefs.h"
#include "statistics.h"
#include "paths_and_components.h"
#include "graph_properties.h"
#include "graph.h"

using namespace std;


namespace CDLib 
{
    string get_graph_details(const graph& g);
    string get_graph_details_components(const graph& g);
    string get_graph_details_efficiency(const graph& g);
    string get_graph_comparisons(const graph& g);
};


#endif	/* GRAPH_SUMMARY_H */

