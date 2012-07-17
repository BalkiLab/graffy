/* 
 * File:   local_community.h
 * Author: prashant
 *
 * Created on June 16, 2012, 12:09 PM
 */

#ifndef LOCAL_COMMUNITY_H
#define	LOCAL_COMMUNITY_H

#include "graph.h"
#include "random.h"
using namespace std;



namespace CDLib
{
    /*
     * Implements Aaron Clauset, "Finding Local Community Structure in Networks" 2005
     * 
     * Input: The graph g
     *        The seed vertex src
     *        The size of the community desired k
     */
    
    void local_community_clauset(const graph& g, id_type src,size_t k, vector< pair<id_type,double> >& output);
    void local_community_clauset_modified(const graph& g, id_type src, size_t k, vector< pair<id_type,double> >& output);
    void LWP_2006(const graph& g, id_type src, vector<id_type>& output);
    void VD_2011(const graph& g, id_type src, node_set& output,id_type k);
    pair<double,node_set> Bagrow_2007(const graph& g, id_type src);
    void My_Algorithm(const graph& g, id_type src, node_set& output);
}


#endif	/* LOCAL_COMMUNITY_H */

