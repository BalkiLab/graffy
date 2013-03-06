/* 
 * File:   centrality.cpp
 * Author: bharath
 * 
 * Created on 5 April, 2012, 2:36 AM
 */

#include "centrality.h"
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
//    Parallel calculation of betweenness centrality. 
//    Parallelism is automatically controlled to number of available CPU's
//    Mutual Exclusion of variable bc is not required as only one addition is performed on that variable.
#ifdef ENABLE_MULTITHREADING
#pragma omp parallel for
#endif
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
    if(!g.is_directed())
#ifdef ENABLE_MULTITHREADING
        #pragma omp parallel for
#endif
        for(id_type i=0;i<bc.size();i++) bc[i] /= 2;
}

void CDLib::betweenness_centralities_normalized(const graph& g, vector<double>& bcn)
{
    betweeness_centralities(g,bcn);
    if (g.get_num_nodes() <= 1)
        return;
    double sum_norm = accumulate(bcn.begin(),bcn.end(),0.0);
    for (id_type i=0; i < bcn.size(); i++)
        bcn[i] /= sum_norm;
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

void CDLib::degree_centralities(const graph& g, vector<id_type>& degrees)
{
//    Return the out-degree centrality of all the nodes in the given graph.
    degrees.clear();
    for(id_type i =0; i < g.get_num_nodes(); i++)
        degrees.push_back(g.get_node_out_degree(i));
}

void CDLib::degree_centralities_normalized(const graph& g, vector<double>& degrees)
{
//    Return the normalized out-degree centrality of all the nodes in the given graph.
    degrees.clear();
    for(id_type i =0; i < g.get_num_nodes(); i++)
        degrees.push_back((double)g.get_node_out_degree(i)/g.get_num_edges());
}

//Overloaded function for a single node
double CDLib::node_clustering_coefficient(const graph&g, id_type node)
{
    if (g.get_node_out_degree(node) == 0)
        return 0;
    if (g.get_node_out_degree(node) == 1)
        return 1;
    
    double edge_count = 0;
    node_set neighbours;
    neighbours.insert(node);
    for(adjacent_edges_iterator aeit = g.out_edges_begin(node);aeit != g.out_edges_end(node);aeit++)
        neighbours.insert(aeit->first);
    for(node_set::iterator nit = neighbours.begin(); nit != neighbours.end(); nit++){
        for(adjacent_edges_iterator aeit = g.out_edges_begin(*nit);aeit != g.out_edges_end(*nit); aeit++){
            if(neighbours.find(aeit->first)!=neighbours.end())
                edge_count++;
        }
    }
//    This below line automatically takes care of directed and undirected graphs.
    edge_count /= ((g.get_node_out_degree(node) + 1) * g.get_node_out_degree(node));
    return edge_count;
}

//Overloaded funtion for all nodes.
void CDLib::node_clustering_coefficient(const graph&g, vector<double>& nodes)
{
    nodes.clear();
    nodes.assign(g.get_num_nodes(),0);
#ifdef ENABLE_MULTITHREADING
        #pragma omp parallel for
#endif
    for(id_type i=0; i<g.get_num_nodes();i++)
        nodes[i] = node_clustering_coefficient(g,i);
}

void CDLib::node_clustering_coefficient_normalized(const graph& g, vector<double>& nodes)
{
    node_clustering_coefficient(g,nodes);
    if (g.get_num_nodes() <= 1)
        return;
    double sum_norm = accumulate(nodes.begin(),nodes.end(),0.0);
    for (id_type i=0; i < nodes.size(); i++)
        nodes[i] /= sum_norm;
}

double CDLib::closeness_centrality_original(const graph& g, id_type node)
{
    if ((node <= g.get_num_nodes()) && (node >= 0)){
        return -1;      // Reporting node out of range error.
    }
    vector<double> distances;
    vector< vector<id_type> > preds;
    single_source_shortest_paths_djikstra(g,node,distances,preds);
    double sum = 0;
    for (id_type i=0; i<distances.size(); i++)
        sum += distances[i];
    sum = g.get_num_nodes()/sum;
    return sum;
}

void closeness_centralities_original(const graph& g, vector<double>& closeness)
{
    closeness.clear();
    closeness.assign(g.get_num_nodes(),0);
    vector< vector<double> > distance_matrix;
    all_pairs_shortest_paths(g,distance_matrix);
    for(id_type i=0; i<g.get_num_nodes();i++){
        double sum = 0;
        for (id_type j=0; j<distance_matrix[i].size(); j++)
                sum += distance_matrix[i][j];
        closeness[i] = g.get_num_nodes()/sum;
    }
}

double CDLib::closeness_centrality(const graph& g, id_type node)
{
//    This approach is explained in the book: Networks, An Introduction by Newman.
//    It takes the Harmonic Mean of the Geodesic Distance.
    if ((node <= g.get_num_nodes()) && (node >= 0)){
        return -1;      // Reporting node out of range error.
    }
    vector<double> distances;
    vector< vector<id_type> > preds;
    single_source_shortest_paths_djikstra(g,node,distances,preds);
    double sum = 0;
    for (id_type i=0; i<distances.size(); i++)
        if (node != i)
            sum += 1/distances[i];
    sum /= g.get_num_nodes() - 1;
    return sum;
}

void CDLib::closeness_centralities(const graph& g, vector<double>& closeness)
{
//    This approach is explained in the book: Networks, An Introduction by Newman.
//    It takes the Harmonic Mean of the Geodesic Distance.
    closeness.clear();
    closeness.assign(g.get_num_nodes(),0);
    vector< vector<double> > distance_matrix;
    all_pairs_shortest_paths(g,distance_matrix);
    for(id_type i=0; i<g.get_num_nodes();i++){
        double sum = 0;
        for (id_type j=0; j<distance_matrix[i].size(); j++)
            if (i != j)
                sum += 1/distance_matrix[i][j];
        closeness[i] = sum/(g.get_num_nodes() - 1);
    }
}

void CDLib::closeness_centralities_normalized(const graph& g, vector<double>& closeness)
{
    closeness_centralities(g, closeness);
    if (g.get_num_nodes() <= 1)
        return;
    double sum_norm = accumulate(closeness.begin(),closeness.end(),0.0);
    for (id_type i=0; i < closeness.size(); i++)
        closeness[i] /= sum_norm;
}

double sum_of_squares (double x, double y) {return (x +(y*y));}

void CDLib::eigenvector_centralities(const graph& g, vector<double>& eigenvector)
{
    eigenvector.clear();
    if (g.get_num_nodes() <= 0)
        return;
    if (g.get_num_nodes() == 1){
        eigenvector.push_back(1);
        return;
    }
    double init, sum_outv, sum_inv;
    sum_inv = sqrt(g.get_num_nodes());
    init = 1/sum_inv;
    sum_outv = 0;
    vector<double> invector(g.get_num_nodes(),init);
    vector<double> outvector;
    id_type iteration_count = 0;        // This is to forcefully terminate the loop depending on # of iteration.
    const id_type converge = 1000;      // This is max. # of iteration before termination of the loop forcefully.
    while ((abs(sum_inv - sum_outv) > 0.1) && (iteration_count < converge))
    {
        iteration_count++;
        outvector.clear();
        outvector.assign(g.get_num_nodes(),0);
#pragma omp parallel for
        for(id_type i=0;i<g.get_num_nodes();i++){
            for(adjacent_edges_iterator aeit = g.out_edges_begin(i);aeit != g.out_edges_end(i);aeit++){
                outvector[i] += aeit->second * invector[aeit->first];
            }
        }
        double norm = accumulate(outvector.begin(),outvector.end(),0.0,sum_of_squares);
        if (norm > 0)
            norm = sqrt(norm);
        sum_inv = 0;            sum_outv = 0;
#pragma omp parallel for
        for(id_type i=0;i<g.get_num_nodes();i++){
            sum_inv += invector[i];
            invector[i] = outvector[i]/norm;
            sum_outv += invector[i];
        }
    }
    eigenvector.assign(invector.begin(),invector.end());
}

void CDLib::eigenvector_centralities_normalized(const graph& g, vector<double>& eigenvector)
{
    eigenvector_centralities(g, eigenvector);
    if (g.get_num_nodes() <= 1)
        return;
    double sum_norm = accumulate(eigenvector.begin(),eigenvector.end(),0.0);
    for (id_type i=0; i < eigenvector.size(); i++)
        eigenvector[i] /= sum_norm;
}