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
    
    bool local_community_clauset(const graph& g, id_type src,size_t k,node_set& output);
    bool local_community_clauset_modified(const graph& g, id_type src, size_t k, node_set& output);
    bool LWP_2006(const graph& g, id_type src, node_set& output);
    void VD_2011(const graph& g, id_type src, node_set& output);
    bool Bagrow_2007(const graph& g, id_type src, node_set& output);
    bool My_Algorithm(const graph& g, id_type src, node_set& output);
    bool CZR(const graph& g, id_type src, node_set& output);
    bool CZR_Beta(const graph& g, id_type src, node_set& output);
}


#endif	/* LOCAL_COMMUNITY_H */

