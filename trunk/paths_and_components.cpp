/* 
 * File:   community_tools.cpp
 * Author: bharath
 * 
 * Created on April 2, 2012, 3:18 PM
 */


#include "paths_and_components.h"

using namespace CDLib;

void CDLib::dfs_visitor(const graph& g, node_set& visited, id_type source)
{
    visited.insert(source);
    for(adjacent_edges_iterator aeit = g.out_edges_begin(source); aeit != g.out_edges_end(source); aeit++)
        if(visited.find(aeit->first) == visited.end()) dfs_visitor(g,visited,aeit->first);
}

void CDLib::bfs_visitor(const graph& g, node_set& visited, id_type source)
{
    queue<id_type> q_bfs;
    q_bfs.push(source);
    while(!q_bfs.empty())
    {
        id_type current  = q_bfs.front();
        q_bfs.pop();
        visited.insert(current);
        for(adjacent_edges_iterator aeit = g.out_edges_begin(current);aeit != g.out_edges_end(current);aeit++)
        {
            if(visited.find(aeit->first) == visited.end())
            {
                visited.insert(aeit->first);
                q_bfs.push(aeit->first);
            }
        }
    }
}

double CDLib::single_source_shortest_paths_bfs(const graph& g,id_type source,vector<double>& distances,vector< vector<id_type> >& preds)
{
    distances.assign(g.get_num_nodes(),numeric_limits<double>::infinity());
    preds.assign(g.get_num_nodes(),vector<id_type>());
    queue<id_type> q_bfs;
    q_bfs.push(source);
    distances[source] = 0;
    id_type last_node = source;
    while(!q_bfs.empty())
    {
        id_type current  = q_bfs.front();
        last_node = current;
        q_bfs.pop();
        for(adjacent_edges_iterator aeit = g.out_edges_begin(current);aeit != g.out_edges_end(current);aeit++)
        {
            if(distances[aeit->first] == numeric_limits<double>::infinity())
            {
                distances[aeit->first] = distances[current]+ 1;
                q_bfs.push(aeit->first);
                
            }
            if(distances[aeit->first] == distances[current]+ 1)
                preds[aeit->first].push_back(current);
        }
    }
    return distances[last_node];
}

bool CDLib::is_path_present(const graph& g, id_type source, id_type dest)
{
    node_set visited;
    bfs_visitor(g,visited,source);
    return (visited.find(dest) != visited.end());
}

id_type CDLib::get_component_around_node_undirected(const graph& g,id_type id, node_set& visited)
{
    if(g.is_directed()) return 0;
    bfs_visitor(g,visited,id);
    return visited.size();
}

id_type CDLib::get_component_around_node_weak(const graph&g, id_type id,node_set& visited)
{
    if(g.is_directed())
    {
        graph g_temp(g);
        g_temp.convert_to_undirected();
        return get_component_around_node_undirected(g_temp,id,visited);
    }
    return get_component_around_node_undirected(g,id,visited);
}

void tarjan_strongconnect(const graph& g,id_type curr,node_set& visited,vector<id_type>& indices,vector<id_type>& lowlinks)
{
    static id_type index = 0;
    indices[curr] = index;
    lowlinks[curr] = index;
    index++;
    visited.insert(curr);
    for(adjacent_edges_iterator aeit = g.out_edges_begin(curr);aeit != g.out_edges_end(curr);aeit++)
    {
        if(!indices[aeit->first])
        {
            tarjan_strongconnect(g,aeit->first,visited,indices,lowlinks);
            lowlinks[curr] = min(lowlinks[curr],lowlinks[aeit->first]);
        }
        else if(visited.find(aeit->first) != visited.end())
            lowlinks[curr] = min(lowlinks[curr],indices[aeit->first]);
    }
    if(indices[curr] != lowlinks[curr])
    {
        visited.clear();
        visited.insert(curr);
    }
}

id_type CDLib::get_component_around_node_strong(const graph& g, id_type id, node_set& visited)
{
    if(g.is_directed())
    {
        vector<id_type> indices(g.get_num_nodes(),0),lowlinks(g.get_num_nodes(),0);
        tarjan_strongconnect(g,id,visited,indices,lowlinks);
        return visited.size();
    }
    return get_component_around_node_undirected(g,id,visited);
}

bool CDLib::is_connected_undirected(const graph& g)
{
    if(g.is_directed()) return false;
    node_set visited;
    return (get_component_around_node_undirected(g,0,visited) == g.get_num_nodes());
}
bool CDLib::is_connected_weakly(const graph& g)
{
    if(!g.is_directed()) return is_connected_undirected(g);
    node_set visited;
    return (get_component_around_node_weak(g,0,visited) == g.get_num_nodes());
}
bool CDLib::is_connected_strongly(const graph& g)
{
     if(!g.is_directed()) return false;
     node_set visited;
     return (get_component_around_node_strong(g,0,visited) == g.get_num_nodes());
}

id_type CDLib::get_largest_connected_component(const graph& g, node_set &members)
{
//Return the size of the largest connected component in the graph. It also 
//populates the set 'members' with the members of the largest connected component.
//Returns strongly connected components for a directed graph.
    vector<node_set> components;
    members.clear();
    if(g.is_directed()){
        get_connected_components_undirected(g,components);
    }
    else{
        get_strongly_connected_components(g,components);
    }
    id_type max_size = 0, max_id = 0;
    for (id_type i=0; i < components.size(); i++){
        if (components[i].size() > max_size){
            max_size = components[i].size();
            max_id = i;
        }
    }
    members.insert(components[max_id].begin(),components[max_id].end());
    return max_size;
}

id_type CDLib::get_connected_components_undirected(const graph& g, vector<node_set>& components)
{
    if(g.is_directed()) return 0;
    components.clear();
    for(id_type i = 0; i<g.get_num_nodes(); i++)
    {
        bool is_found = false;
        for(id_type j=0;j<components.size();j++)
        {
            //is_found = is_found || components[j].find(i) != components[j].end();
            if(components[j].find(i) != components[j].end())
            {
                is_found = true;
                break;
            }
        } 
        if(!is_found)
        {
            components.push_back(node_set());
            get_component_around_node_undirected(g,i,components[components.size()-1]);
        }
    }
    return components.size();
}

id_type CDLib::get_weakly_connected_components(const graph& g, vector<node_set>& components)
{
    if(g.is_directed())
    {
        graph g_temp(g);
        g_temp.convert_to_undirected();
        return get_connected_components_undirected(g_temp,components);
    }
    else return get_connected_components_undirected(g,components);
}

id_type CDLib::get_strongly_connected_components(const graph& g, vector<node_set>& components)
{
    if(g.is_directed())
    {
        for(id_type i = 0; i<g.get_num_nodes(); i++)
        {
            bool is_found = false;
            for(id_type j=0;j<components.size();j++) is_found = is_found || (components[j].find(i) != components[j].end());
            if(!is_found)
            {
                components.push_back(node_set());
                get_component_around_node_strong(g,i,components.back());
            }
        }
         return components.size();
    }
    else return get_connected_components_undirected(g,components);
}

bool visit_topsort(const graph& g,id_type id,vector<id_type>& ordering,node_set& visited)
{
    visited.insert(id);
    for(adjacent_edges_iterator aeit = g.in_edges_begin(id);aeit != g.in_edges_end(id);aeit++)
    {
        if(visited.find(aeit->first) != visited.end()) return false;
        return visit_topsort(g,aeit->first,ordering,visited);
    }
    ordering.push_back(id);
    return true;
}

bool CDLib::get_topological_ordering(const graph& g,vector<id_type>& ordering)
{
    node_set visited;
    for(id_type i=0; i< g.get_num_nodes();i++)
    {
        if(!g.get_node_out_degree(i) && !visit_topsort(g,i,ordering,visited)) return false; 
    }
    return true;
}

double CDLib::single_source_shortest_paths_djikstra(const graph&g,id_type source,vector<double>& distances,vector< vector<id_type> >& preds)
{
    if(!g.is_weighted()) single_source_shortest_paths_bfs(g,source,distances,preds);
    else if(!has_negative_edge_weights(g))
    {
        distances.assign(g.get_num_nodes(),numeric_limits<double>::infinity());
        preds.assign(g.get_num_nodes(),vector<id_type>());
        distances[source] = 0;
        binary_heap p_queue(distances,false);
        id_type last_node = source;
        while(!p_queue.empty())
        {
            id_type curr = p_queue.top().first;
            last_node = curr;
            p_queue.pop();
            if(distances[curr] == numeric_limits<double>::infinity()) break;
            for(adjacent_edges_iterator aeit = g.out_edges_begin(curr); aeit != g.out_edges_end(curr);aeit++)
            {
                double alt = distances[curr] + aeit->second;
                if(alt <= distances[aeit->first])
                {
                    distances[aeit->first] = alt;
                    p_queue.update_key(make_pair(aeit->first,alt));
                    preds[aeit->first].push_back(curr);
                }
            }
        }
        return distances[last_node];
    }
    return 0;
}

double CDLib::diameter(const graph& g)
{
    double max =0;
    for(id_type i=0;i<g.get_num_nodes();i++)
    {
        vector<double> distances;
        vector< vector<id_type> > preds;
        double test_diameter;
        if(!g.is_weighted()) test_diameter = single_source_shortest_paths_bfs(g,i,distances,preds);
        else test_diameter = single_source_shortest_paths_djikstra(g,i,distances,preds);
        if(max < test_diameter) max = test_diameter;
    }
    return max;
}

void CDLib::all_pairs_shortest_paths(const graph& g, vector< vector<double> >& path_matrix)
{
//    The node_control value is set to control the memory overload due to large number of nodes. This value can be
//    set depending on the available memory of the system. The memory calculation is:
//    memory = node_control x node_control x 8 x 2 + some small constant
    id_type node_control = 10000;
    if (g.get_num_nodes() < 1 )
        return;
    else if (g.get_num_nodes() == 1){
        path_matrix.clear();
        path_matrix.assign(g.get_num_nodes(),vector<double>(g.get_num_nodes(),0));
        return;
    }
    else if (g.get_num_nodes() < node_control){
        if (g.get_density() < (log(2)/(log(g.get_num_nodes()))))
            all_pairs_shortest_paths_djikshtra(g,path_matrix);
        else
            all_pairs_shortest_paths_floyd_warshal(g,path_matrix);
    }
    else{
        all_pairs_shortest_paths_djikshtra(g,path_matrix);
    }
}

void CDLib::all_pairs_shortest_paths_djikshtra(const graph& g, vector< vector<double> >& path_matrix)
{
    path_matrix.assign(g.get_num_nodes(),vector<double>());
#pragma omp parallel for     
    for(id_type i=0;i<g.get_num_nodes();i++)
    {
        vector< vector<id_type> > preds;
        if(!g.is_weighted()) single_source_shortest_paths_bfs(g,i,path_matrix[i],preds);
        else single_source_shortest_paths_djikstra(g,i,path_matrix[i],preds);
    }
}

void CDLib::all_pairs_shortest_paths_floyd_warshal(const graph& g, vector< vector<double> >& path_matrix)
{
//    Should be used only when space is not a constraint.
    path_matrix.clear();
    vector< vector<double> > prev (g.get_num_nodes(),vector<double>(g.get_num_nodes(),numeric_limits<double>::infinity()));
    vector< vector<double> > next (g.get_num_nodes(),vector<double>(g.get_num_nodes(),numeric_limits<double>::infinity()));

    
    for(id_type i=0;i<g.get_num_nodes();i++){
        for(adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i);aeit++){
            prev[i][aeit->first] = aeit->second;
        }
        prev[i][i] = 0;
    }
    for(id_type k=0;k<g.get_num_nodes();k++){
#pragma omp parallel for        
        for(id_type i=0;i<g.get_num_nodes();i++){
            for(id_type j=0;j<g.get_num_nodes();j++){
                next[i][j] = ((prev[i][j] < (prev[i][k] + prev[k][j])) ? prev[i][j] : (prev[i][k] + prev[k][j]));
            }
        }
        prev.clear();
        prev.assign(next.begin(),next.end());
        next.clear();
        if (k < (g.get_num_nodes() - 1))
            next.assign(g.get_num_nodes(),vector<double>(g.get_num_nodes(),numeric_limits<double>::infinity()));
    }
    path_matrix.assign(prev.begin(),prev.end());
}

void CDLib::single_source_shortest_paths_djikstra_with_paths(const graph&g,id_type source,vector<double>& distances,vector< vector<id_type> >& paths)
{
    vector< vector<id_type> > preds;
    single_source_shortest_paths_djikstra(g,source,distances,preds);
    paths.assign(g.get_num_nodes(),vector<id_type>());
    for (id_type i=0; i<g.get_num_nodes(); i++)
    {
        if (preds[i].size() !=0)
        {
            vector<id_type> reverse;
            id_type curr = i;
            reverse.push_back(curr);
            while (preds[curr][0] != source)
            {
                curr = preds[curr][0];
                reverse.push_back(curr);
            }
            paths[i].push_back(source);
            vector<id_type>::reverse_iterator rit;
            for (rit=reverse.rbegin() ; rit < reverse.rend(); ++rit )
                paths[i].push_back(*rit);
            
        }
    }
}

bool CDLib::has_negative_edge_weights(const graph& g)
{
    for(id_type i=0; i<g.get_num_nodes();i++)
        for(adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i);aeit++)
            if(aeit->second < 0) return true;
    return false;
}

void do_count_paths(const graph& g, id_type source, id_type dest, unordered_set<id_type>& visited,vector<id_type>& paths)
{
	for(adjacent_edges_iterator aeit = g.out_edges_begin(source);aeit != g.out_edges_end(source); aeit++)
	{
		if(aeit->first == dest) paths[visited.size()]++;
		else if(visited.find(aeit->first) == visited.end())
		{
			unordered_set<id_type> visited2(visited);
			visited2.insert(aeit->first);
			do_count_paths(g,aeit->first,dest,visited2,paths);
		}
	}
}

void CDLib::get_all_paths(const graph& g,id_type source, id_type dest,vector<id_type>& paths)
{
    paths.clear();
    paths.assign(g.get_num_nodes(),0);
    unordered_set<id_type> visited;
    visited.insert(source);
    do_count_paths(g,source,dest,visited,paths);
}

void CDLib::all_path_lenth_Monte_Carlo(const graph& g, vector< vector<double> >& paths, long monte_c)
{
    // Returns the path length in terms of hop distance using random walk from source-destination.
    /* less than 50% not reachable => avg. path length of that else => not reachable
    * monte_c tell us the number of iteration over which the the mean is taken. */
    vector< vector<id_type> > mc_monte_times (g.get_num_nodes(),vector<id_type>(g.get_num_nodes(),0));
    paths.assign(g.get_num_nodes(),vector<double>(g.get_num_nodes(),0.0));
    for (long s=0; s<monte_c; s++)
    {
        #pragma omp parallel for
        for(id_type i=0; i<g.get_num_nodes();i++)
        {
            id_type x = i,last = i, length = 0;
            paths[i][x] += length;  mc_monte_times[i][x]++;
            graph g2(g);
            g2.convert_to_directed();
            if (g2.get_node_out_degree(x) >= 1)
            {
                length++;
                RandomGenerator<id_type> gen(0,(g2.get_node_out_degree(x)-1));
                id_type assign_index = gen.next();
                id_type new_node,count=0;
                node_set adj;
                for(adjacent_edges_iterator aeit = g2.out_edges_begin(x);aeit != g2.out_edges_end(x); aeit++)
                    adj.insert(aeit->first);
                for (node_set::iterator sit = adj.begin(); sit != adj.end(); sit++, count++)
                {
                    g2.remove_edge(*sit,x);
                    paths[i][*sit] += length;  mc_monte_times[i][*sit]++;
                    if (count == assign_index)
                        new_node = *sit;
                    for (node_set::iterator sjt = adj.begin(); sjt != adj.end(); sjt++)
                    {
                        g2.remove_edge(*sit,*sjt);
                        g2.remove_edge(*sjt,*sit);
                    }
                }
                x = new_node;
            }
            while(g2.get_node_out_degree(x) > 1)
            {
                g2.remove_edge(x,last);
                g2.remove_edge(last,x);
                RandomGenerator<id_type> gen(0,(g2.get_node_out_degree(x)-1));
                id_type assign_index = gen.next();
                id_type new_node,count = 0;
                node_set adj;
                for(adjacent_edges_iterator aeit = g2.out_edges_begin(x);aeit != g2.out_edges_end(x); aeit++)
                    adj.insert(aeit->first);
                for (node_set::iterator sit = adj.begin(); sit != adj.end(); sit++, count++)
                {
                    g2.remove_edge(*sit,x);
                    if (count == assign_index)
                        new_node = *sit;
                }
                last = x;
                x= new_node;
                paths[i][x] += ++length;  mc_monte_times[i][x]++;
            }
            
        }
    }
    unsigned long thres = monte_c * 0.1;
    for(id_type i=0; i<g.get_num_nodes();i++)
    {
        for(id_type j=0; j<g.get_num_nodes();j++)
        {     
            if (mc_monte_times[i][j] >= thres)
            {
                paths[i][j] /= mc_monte_times[i][j];
            }
            else
                paths[i][j] = numeric_limits<double>::infinity();
        }
    }
}

double CDLib::blocking_probability(id_type number_of_nodes, id_type degree, id_type visited)
{
    if ((visited > 1) && (degree == 1))
        return 1;
    else 
        if (visited <= (degree + 1))
        return 0;
    else
    {
        id_type con0 = degree + 1;
        id_type c1 = number_of_nodes - con0;
        id_type c2 = number_of_nodes - 2;
        id_type c4 = visited - 2;
        id_type c3 = visited - con0;
        
        double con4 = c4*log(c4);
        double con1 = c1*log(c1);
        double con3 = c3*log(c3);
        double con2 = c2*log(c2);
        
        double const1 = sqrt((c4*c1)/(c3*c2));
        double const2 = exp(con4 + con1 - con3 - con2);
        return (const1 * const2);
    }
}

void CDLib::alternate_path_length_destabilization(graph&g,id_type source,vector<double>& alternate_distances)
{
    /* This function returns alternate distance (not geodesic) from source to all other nodes.
     * The alternate distance is obtained as an average geodesic distance between a pair of nodes 
     * when each edge in the geodesic path is removed one by one barring from the edge which 
     * would disconnect the graph.*/
    vector<double> dist;
    vector< vector<id_type> > paths;
    single_source_shortest_paths_djikstra_with_paths(g,source,dist,paths);
    alternate_distances.assign(g.get_num_nodes(),0);
    for(id_type i=0; i<g.get_num_nodes();i++)
    {
        if ((dist[i] <=2) || (dist[i] == numeric_limits<double>::infinity()))
            alternate_distances[i] = dist[i];
        else
        {
            long inf = 0;
            for (unsigned long y=0; y<(paths[i].size()-1); y++)
            {
                g.remove_edge(paths[i][y],paths[i][y+1]);
                vector<double> new_dist;
                vector< vector<id_type> > new_preds;
                single_source_shortest_paths_djikstra(g,source,new_dist,new_preds);
                if (new_dist[i] == numeric_limits<double>::infinity())
                    inf++;
                else
                    alternate_distances[i] += new_dist[i];
                g.add_edge(paths[i][y],paths[i][y+1],1);
            }
            id_type operations = paths[i].size() - 1 - inf;
            if (operations == 0)
                alternate_distances[i] = dist[i];
            else
                alternate_distances[i] /= operations;
        }
    }
}

double CDLib::efficiency_sw_global(const graph& g, bool type)
{
    /* Returning Global Efficiency of a Small World Network according to 2001 paper 
     * The bool type represents the type of distance calculation of efficiency. 
     * type = 0 implies distance calculation based on dijkshtra's algorithm. 
     * type = 1 implies distance calculation based on alternate_path_length_destabilization()
     *  function in this library. */
    if (g.get_num_nodes() > 1)
    {
        double efficiency = 0;
        double ideal = 1 / g.minimum_weight();
        for (unsigned long i = 0 ; i < g.get_num_nodes(); i++)
        {
            vector< vector<id_type> > preds;
            vector<double> distance;
            if (type){
                graph gtemp(g);
                alternate_path_length_destabilization(gtemp,i,distance);
            }
            else
                single_source_shortest_paths_djikstra(g,i,distance,preds);
            for (unsigned long j = 0 ; j < g.get_num_nodes(); j++)
                if (i != j)
                    efficiency = efficiency + (1/distance[j]);
        }
//        cout << "\nIdeal :" << ideal << "\tEfficiency :" << efficiency << endl;
        efficiency = efficiency / (g.get_num_nodes() * (g.get_num_nodes() - 1));
        efficiency /= ideal;
        return efficiency;
    }
    else if (g.get_num_nodes() == 1)
            return 1;
    else
            return 0;
}

double CDLib::efficiency_sw_global_monte_carlo(graph& g)
{
    /* Returns the small world global efficiency of the graph using loop avoiding 
     * random walk and averaging out using monte carlo simulation as the distance 
     * measure in the efficiency formula. The monte carlo simulation is parallelized
     * using open MPI #prgma. */
    vector< vector<double> > paths;
    double efficiency = 0;
    all_path_lenth_Monte_Carlo(g, paths, 10000);
    for (unsigned long i = 0 ; i < g.get_num_nodes(); i++)
        for (unsigned long j = 0 ; j < g.get_num_nodes(); j++)
            if (i != j)
                    efficiency += 1 / paths[i][j];
    efficiency /= g.get_num_nodes() * (g.get_num_nodes() - 1);
    return efficiency;
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

double CDLib::path_entropy(graph& g)
{
    // Centrality Entropy of the Graph based on shortest path connectivity, as in Borgatti paper
    double entropy = 0;
    vector<id_type> each_source (g.get_num_nodes(),0);
    id_type all_sources = 0;
    for (id_type i = 0; i < g.get_num_nodes(); i++){
        vector<double> distances;
        vector< vector<id_type> > preds;
        single_source_shortest_paths_bfs(g,i,distances,preds);
        for (id_type j = 0; j < g.get_num_nodes(); j++){
            if ((j!=i) && (distances[j] != numeric_limits<double>::infinity()) && (distances[j] != 0)){
                each_source[i]++;
                all_sources++;
            }
        }
    }
    if (all_sources > 0){
        for (id_type i = 0; i < g.get_num_nodes(); i++){
            double prob = (double)each_source[i]/all_sources;
            if (prob != 0){ entropy += prob * (log(prob)/log(2)); }
        }
        entropy *= -1;
    }
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

/*degree_dist_controlling_parameter controls degree distribution between power law and exponential degree distribution as it varies between 0 to 1.
  it is undirected graph model.*/

double energy(graph& g,double degree_dist_controlling_parameter,id_type nC2,double Dlinear);

void CDLib::generate_ferrer_i_cancho_model(graph& g,size_t num_nodes, size_t max_failure_allowed,double degree_dist_controlling_parameter,double probability_to_alter_edge,double initial_probability_of_edge)
{
    if(degree_dist_controlling_parameter<=1 && degree_dist_controlling_parameter>=0 && initial_probability_of_edge<=1 && initial_probability_of_edge>=0 && probability_to_alter_edge<=1 && probability_to_alter_edge>=0)
    {
        double Dlinear = ((num_nodes + 1) / 3.0);
        id_type nC2 = ((num_nodes * (num_nodes - 1)) / 2);
        
        UniformRandomGenerator<double> rand;

        graph g_t(0, 0);
        init_empty_graph(g,num_nodes);
        init_empty_graph(g_t,num_nodes);
      
        for(id_type i=0;i<num_nodes;i++)
        {
            for(id_type j=i+1;j<num_nodes;j++)
            {
                double R2=rand.next(1);
                if(R2<initial_probability_of_edge && i!=j)
                {
                   g.add_edge(i,j,1);
                } 
            }
        }

        id_type failure = 0;

        double Energy_old = energy(g, degree_dist_controlling_parameter, nC2, Dlinear);
        
        while (failure < max_failure_allowed) 
        {
            g_t=g;

            for (id_type i = 0; i < num_nodes; i++) 
            {
                for (id_type j = i+1 ; j < num_nodes; j++) 
                {
                    double R1=rand.next(1);
                    if(R1<probability_to_alter_edge)
                    {
                        if(g_t.get_edge_weight(i,j)==0)
                            g_t.add_edge(i,j,1);
                        else
                            g_t.remove_edge(i,j);
                    }
                }
            }

            double Energy_new = energy(g_t, degree_dist_controlling_parameter, nC2, Dlinear);

            if (Energy_new < Energy_old) 
            {
                g = g_t;
                Energy_old = Energy_new;
                failure = 0;
            } else
                failure++;
        }
    }
    else
        cout<<"\n last three parameter value should be between 0 and 1";
    
}

double energy(graph& g,double degree_dist_controlling_parameter,id_type nC2,double Dlinear)
{
    id_type num_edges = g.get_num_edges();
    id_type num_nodes = g.get_num_nodes();
    
    vector<double> distance;
    vector< vector <id_type> > seq;
    
    id_type total_of_min_distance=0;
    
    for(id_type i=0;i<num_nodes;i++)
    {
        single_source_shortest_paths_bfs(g,i,distance,seq);
        for(id_type j=0;j<num_nodes;j++)
        {
            total_of_min_distance = total_of_min_distance + distance[j];
        }
    }
    double normalized_num_of_link = (double)num_edges/nC2;
    double avg_min_distance = (double)total_of_min_distance/(nC2*2);
    
    double normalized_distance = avg_min_distance/Dlinear;
    double Energy = degree_dist_controlling_parameter*normalized_distance + (1-degree_dist_controlling_parameter)*normalized_num_of_link;
    
    return Energy;
}

id_type CDLib::hop_distance_matrix(const graph& g, vector< vector<id_type> > & path_matrix)
{
    id_type max_dist = 0;
    path_matrix.assign(g.get_num_nodes(),vector<id_type>(g.get_num_nodes(),numeric_limits<id_type>::max()));
    for(id_type i=0;i<g.get_num_nodes();i++)
    {
        path_matrix[i][i] = 0;
        vector<double> dist;
        vector< vector<id_type> > preds;
        single_source_shortest_paths_bfs(g,i,dist,preds);
        for(id_type j=0;j<dist.size();j++)
        {
            path_matrix[i][j] = dist[j];
            if(dist[j]>max_dist) max_dist = dist[j];
        }
    }
    return max_dist;
}