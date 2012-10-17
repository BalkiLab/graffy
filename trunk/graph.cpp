/* 
* File:   graph.cpp
* Author: bharath
* 
* Created on April 1, 2012, 12:00 PM
*/

#include "graph.h"
using namespace CDLib;

graph::graph(): b_directed(false),b_weighted(false), blm_labels(),dam_backend(){}
graph::graph(bool directed, bool weighted): b_directed(directed),b_weighted(weighted), blm_labels(),dam_backend(){}

bool graph::is_directed() const { return b_directed;}
bool graph::is_weighted() const { return b_weighted; }
       
id_type graph::get_num_nodes() const { return dam_backend.num_nodes();}
id_type graph::get_num_edges() const
{
   if(is_directed()) return dam_backend.num_edges();
   else return (dam_backend.num_edges()+dam_backend.num_self_edges())/2;
}

id_type graph::get_num_self_edges() const { return dam_backend.num_self_edges();}
        
wt_t graph::get_total_weight() const
{
   if(is_directed()) return dam_backend.total_weight();
   else return (dam_backend.total_weight()+dam_backend.self_edges_weight())/2;
}
wt_t graph::get_self_edges_weight() const { return dam_backend.self_edges_weight();}

wt_t graph::get_density() const 
{ 
    id_type n = get_num_nodes();
    double denom = n*((double)n + (double)(get_num_self_edges()? 0:-1));
    return (is_directed()? 1:2)*get_num_edges()/denom;
}

string graph::get_node_label(id_type id) const { return blm_labels.get_label(id); }
id_type graph::get_node_id(const string& label) const { return blm_labels.get_id(label); }

id_type graph::get_node_in_degree(id_type id) const { return dam_backend.in_degree(id);}
id_type graph::get_node_in_degree(const string& label) const { return dam_backend.in_degree(blm_labels.get_id(label));}   

id_type graph::get_node_out_degree(id_type id) const { return dam_backend.out_degree(id);}
id_type graph::get_node_out_degree(const string& label) const { return dam_backend.out_degree(blm_labels.get_id(label));}    

wt_t graph::get_node_in_weight(id_type id) const { return dam_backend.in_degree_wt(id);}
wt_t graph::get_node_in_weight(const string& label) const { return dam_backend.in_degree_wt(blm_labels.get_id(label));}   

wt_t graph::get_node_out_weight(id_type id) const { return dam_backend.out_degree_wt(id);}
wt_t graph::get_node_out_weight(const string& label) const { return dam_backend.out_degree_wt(blm_labels.get_id(label));}   

wt_t graph::get_edge_weight(id_type from_id, id_type to_id) const { return dam_backend.edge_weight(from_id,to_id); }
wt_t graph::get_edge_weight(const string& from_label, const string& to_label) const { return dam_backend.edge_weight(blm_labels.get_id(from_label),blm_labels.get_id(to_label)); }

node_label_iterator graph::node_labels_begin() const { return blm_labels.begin();}
node_label_iterator graph::node_labels_end() const { return blm_labels.end();}

adjacent_edges_iterator graph::in_edges_begin(id_type id) const { return dam_backend.in_edges_begin(id); }
adjacent_edges_iterator graph::in_edges_begin(const string& label) const { return dam_backend.in_edges_begin(blm_labels.get_id(label)); }

adjacent_edges_iterator graph::in_edges_end(id_type id) const { return dam_backend.in_edges_end(id); }
adjacent_edges_iterator graph::in_edges_end(const string& label) const { return dam_backend.in_edges_end(blm_labels.get_id(label)); }

adjacent_edges_iterator graph::out_edges_begin(id_type id) const { return dam_backend.out_edges_begin(id); }
adjacent_edges_iterator graph::out_edges_begin(const string& label) const { return dam_backend.out_edges_begin(blm_labels.get_id(label)); }

adjacent_edges_iterator graph::out_edges_end(id_type id) const { return dam_backend.out_edges_end(id); }
adjacent_edges_iterator graph::out_edges_end(const string& label) const { return dam_backend.out_edges_end(blm_labels.get_id(label)); }

id_type graph::add_node(const string& label)
{
   if(blm_labels.insert(label)) return dam_backend.insert_node();
   return get_num_nodes()-1;
}

id_type graph::add_node()
{
   ostringstream oss;
   oss << get_num_nodes();
   return add_node(oss.str());
}

bool graph::add_edge(id_type from_id, id_type to_id, wt_t weight) 
{
    double weight2 = weight;
    if(!is_weighted() && weight) weight2 = 1;
    if(!is_directed())dam_backend.insert_edge(to_id,from_id,weight2);
    return dam_backend.insert_edge(from_id,to_id,weight2);
}

void graph::add_self_edges(double weight)
{
    if(weight)
        for (id_type i=0; i<get_num_nodes(); i++)
                add_edge(i,i,weight);
}

void graph::remove_self_edges()
{
    for (id_type i=0; i<get_num_nodes(); i++)
            remove_edge(i,i);
}

wt_t graph::add_edge(const string& from_label, const string& to_label,wt_t weight) 
{ 
    return add_edge(blm_labels.get_id(from_label),blm_labels.get_id(to_label),weight);
}

bool graph::remove_node(id_type id) 
{ 
    if(blm_labels.erase(id)) return dam_backend.delete_node(id);
    return false;
}
bool graph::remove_node(const string& label) { return remove_node(get_node_id(label));}

id_type graph::remove_nodes(const set<id_type>& nodes)
{
//  Returns the number of nodes removed.    
    set<string> to_remove;
    for(set<id_type>::iterator it = nodes.begin(); it != nodes.end(); it++)
        to_remove.insert(get_node_label(*it));
    return remove_nodes(to_remove);
}

id_type graph::remove_nodes(const set<string>& nodes)
{
//  Returns the number of nodes removed.    
    id_type number_of_removed = 0;
    for (set<string>::iterator it = nodes.begin(); it != nodes.end(); it++){
        if (remove_node(*it)){ number_of_removed++; }
    }
    return number_of_removed;
}

id_type graph::remove_isolates()
{
// Returns the number of isolates nodes removed.
    id_type number_of_removed = 0;
    set<string> to_remove;
    for(id_type i=0;i<get_num_nodes();i++){
        if ((get_node_in_degree(i) == 0) && (get_node_out_degree(i) == 0)){
            to_remove.insert(get_node_label(i));
        }
    }
    return remove_nodes(to_remove);
}

bool graph::remove_edge(id_type from_id, id_type to_id) 
{ 
    if(!is_directed()) return dam_backend.delete_edge(from_id,to_id) && dam_backend.delete_edge(to_id,from_id);
    else return dam_backend.delete_edge(from_id,to_id);
}
wt_t graph::remove_edge(const string& from_label, const string& to_label) { return remove_edge(blm_labels.get_id(from_label),blm_labels.get_id(to_label)); }

bool graph::set_edge_weight(id_type from_id, id_type to_id, wt_t weight) 
{
    if(!is_weighted()) return dam_backend.set_edge_wt(from_id,to_id,1);
    return dam_backend.set_edge_wt(from_id,to_id,weight);
}
wt_t graph::set_edge_weight(const string& from_label, const string& to_label,wt_t weight) { return set_edge_weight(blm_labels.get_id(from_label),blm_labels.get_id(to_label),weight);}

bool graph::remove_all_edges()
{
    return dam_backend.delete_all_edges();
}

bool graph::clear()
{
    return dam_backend.clear() && blm_labels.clear();
}

bool graph::convert_to_unweighted(double threshold)
{
    if(!is_weighted()) return false;
    for(id_type i=0;i<get_num_nodes();i++)
        for(adjacent_edges_iterator aeit = out_edges_begin(i);aeit != out_edges_end(i);aeit++)
            set_edge_weight(i,aeit->first,((aeit->second > threshold)?1:0));
    return b_weighted = true;
}
bool graph::convert_to_undirected()
{
    if(!is_directed()) return false;
    for(id_type i=0;i<get_num_nodes();i++)
        for(adjacent_edges_iterator aeit = out_edges_begin(i);aeit != out_edges_end(i);aeit++)
            add_edge(aeit->first,i,aeit->second);
    return b_directed = true;
}

bool graph::convert_to_directed()
{
    return b_directed = 1;
}

bool graph::convert_to_weighted()
{
    return b_weighted = 1;
}


double graph::extreme_weight(bool max) const
{
    if(!is_weighted()) return 1;
    double minimum = numeric_limits<double>::infinity();
    double maximum = -1 * numeric_limits<double>::infinity();
    for(id_type i=0;i<get_num_nodes();i++){
        for(adjacent_edges_iterator aeit = out_edges_begin(i);aeit != out_edges_end(i);aeit++){
            if (aeit->second < minimum)
                minimum = aeit->second;
            if (aeit->second > maximum)
                maximum = aeit->second;
        }
    }
    if (max)
        return maximum;
    else
        return minimum;
}

double graph::minimum_weight() const 
{
    return extreme_weight(0);
}

double graph::maximum_weight() const
{
    return extreme_weight(1);
}
