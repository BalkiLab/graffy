/* 
 * File:   centrality.cpp
 * Author: bharath
 * 
 * Created on 5 April, 2012, 2:36 AM
 */

#include "centrality.h"
#include "disjoint_set.h"
using namespace CDLib;

void do_bfs(const graph& g,id_type i,vector< vector<id_type> >& preds,vector<id_type>& paths, vector<double>& dist, stack<id_type>& s_dist)
{
    queue<id_type> q_bfs;
    q_bfs.push(i);
    while(!q_bfs.empty())
    {
        id_type curr = q_bfs.front();
        q_bfs.pop();
        s_dist.push(curr);
        for(adjacent_edges_iterator aeit = g.out_edges_begin(curr);aeit != g.out_edges_end(curr);aeit++)
        {
            if(dist[aeit->first] == numeric_limits<double>::infinity())
            {
                q_bfs.push(aeit->first);
                dist[aeit->first] = dist[curr] + 1;
            }
            if(dist[aeit->first] == dist[curr]+1)
            {
                paths[aeit->first]+= paths[curr];
                preds[aeit->first].push_back(curr);
            }
        }
    }
}

void do_djikstra(const graph& g,id_type i,vector< vector<id_type> >& preds,vector<id_type>& paths, vector<double>& dist, stack<id_type>& s_dist)
{
    binary_heap p_queue(dist,false);
    while(!p_queue.empty())
    {
        id_type curr = p_queue.top().first;
        p_queue.pop();
        if(dist[curr] == numeric_limits<double>::infinity()) break;
        for(adjacent_edges_iterator aeit = g.out_edges_begin(curr); aeit != g.out_edges_end(curr);aeit++)
        {
            double alt = dist[curr] + aeit->second;
            if(dist[aeit->first] == numeric_limits<double>::infinity())
            {
                dist[aeit->first] = alt;
                p_queue.update_key(make_pair(aeit->first,alt));
            }
            if(dist[aeit->first] == alt)
            {
                preds[aeit->first].push_back(curr);
                paths[aeit->first]+= paths[curr];
            }
        }
    }
}

void CDLib::betweeness_centralities(const graph& g, vector<double>& bc)
{
    bc.clear();
    bc.assign(g.get_num_nodes(),0);
    for(id_type i=0; i<g.get_num_nodes();i++)
    {
        vector< vector<id_type> > preds(g.get_num_nodes(),vector<id_type>());
        vector<double> dist(g.get_num_nodes(),numeric_limits<double>::infinity());
        vector<id_type> paths(g.get_num_nodes(),0);
        stack<id_type> s_dist;
        dist[i] = 0;
        paths[i] = 1;
        if(!g.is_weighted()) do_bfs(g,i,preds,paths,dist,s_dist);
        else do_djikstra(g,i,preds,paths,dist,s_dist);
        vector<double> deps(g.get_num_nodes(),0);
        while(!s_dist.empty())
        {
            id_type curr = s_dist.top();
            s_dist.pop();
            for(id_type j=0;j<preds[curr].size();j++)
                deps[preds[curr][j]] += (((double)paths[preds[curr][j]]/(double)paths[curr])*(1+deps[curr]));
            if(curr != i) bc[curr]+=deps[curr];
        }
    }
    if(!g.is_directed())for(id_type i=0;i<bc.size();i++) bc[i] /= 2;
}

double CDLib::edge_clustering_coefficient(const graph&g,id_type from_id, id_type to_id)
{
    if(!g.get_edge_weight(from_id,to_id)) return 0;
    double denom = min(g.get_node_in_weight(from_id)-1,g.get_node_out_weight(to_id)-1);
    if(denom == 0) return numeric_limits<double>::infinity();
    if(!g.is_directed()) denom /=2;
    set<id_type> from_neighs,to_neighs;
    for(adjacent_edges_iterator aeit = g.out_edges_begin(to_id);aeit != g.out_edges_end(to_id);aeit++) from_neighs.insert(aeit->first);
    for(adjacent_edges_iterator aeit = g.in_edges_begin(from_id);aeit != g.in_edges_end(from_id);aeit++) to_neighs.insert(aeit->first);
    vector<id_type> common_neighbors(g.get_node_out_degree(to_id)+g.get_node_in_degree(from_id),0);
    vector<id_type>::iterator common_it = set_intersection(from_neighs.begin(),from_neighs.end(),to_neighs.begin(),to_neighs.end(),common_neighbors.begin());
    double numer = (double)((common_it-common_neighbors.begin()) + 1);
    if(g.is_weighted()) 
    {
        double in_sum = 0,out_sum = 0;
        for(vector<id_type>::iterator it = common_neighbors.begin() ;it != common_it ;it++ )
        {
            out_sum+=g.get_edge_weight(to_id,*it);
            in_sum+=g.get_edge_weight(*it,from_id);
        }
        numer = in_sum*out_sum;
        return numer/denom;
    }
    return (double)numer/denom;
}

/*This gives the #nodes of degree indicated as the index of sequence variable in an undirected graph*/
void CDLib::degree_sequence(const graph& g, vector<id_type>& sequence)
{
    sequence.clear();
    if (g.is_directed())
    {
        id_type max_degree = 0;
        for(id_type i =0; i < g.get_num_nodes(); i++)
            if (max_degree < g.get_node_out_degree(i))
                max_degree = g.get_node_out_degree(i);
        sequence.assign((max_degree + 1),0);
        for(id_type i =0; i < g.get_num_nodes(); i++)
            sequence[g.get_node_out_degree(i)]++;
    }
}