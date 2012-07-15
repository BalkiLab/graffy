/* 
* File:   double_adjacency_map.cpp
* Author: bharath
* 
* Created on April 1, 2012, 1:50 AM
*/

#include "double_adjacency_map.h"
#include "statistics.h"
using namespace CDLib;

double_adjacency_map::double_adjacency_map() : st_num_edges(0), st_num_self_edges(0),wt_total_wt(0),wt_self_edge_wt(0),am_in_edges(),am_out_edges(),vw_in_degree(),vw_out_degree(){}

id_type double_adjacency_map::num_nodes() const {return am_in_edges.size(); }

id_type double_adjacency_map::num_edges() const {return st_num_edges; }

id_type double_adjacency_map::num_self_edges() const { return st_num_self_edges; }

wt_t double_adjacency_map::total_weight() const { return wt_total_wt; }

wt_t double_adjacency_map::self_edges_weight() const { return wt_self_edge_wt; }

bool double_adjacency_map::is_valid_node(id_type id) const { return (id>=0 && id<am_in_edges.size()); }

id_type double_adjacency_map::in_degree(id_type id) const 
{
  if(is_valid_node(id)) return am_in_edges[id].size(); 
  return 0;
}

id_type double_adjacency_map::out_degree(id_type id) const 
{
  if(is_valid_node(id)) return am_out_edges[id].size(); 
  return 0;
} 

wt_t double_adjacency_map::in_degree_wt(id_type id) const 
{ 
  if(is_valid_node(id)) return vw_in_degree[id]; 
  return 0;
}

wt_t double_adjacency_map::out_degree_wt(id_type id) const 
{ 
  if(is_valid_node(id)) return vw_out_degree[id]; 
  return 0;
}   

adjacent_edges_iterator double_adjacency_map::in_edges_begin(id_type id) const
{
  if(is_valid_node(id)) return am_in_edges[id].begin();
  return am_in_edges[am_in_edges.size()-1].end();
}

adjacent_edges_iterator double_adjacency_map::in_edges_end(id_type id) const
{
  if(is_valid_node(id)) return am_in_edges[id].end();
  return am_in_edges[am_in_edges.size()-1].end();
}

adjacent_edges_iterator double_adjacency_map::out_edges_begin(id_type id) const
{
  if(is_valid_node(id)) return am_out_edges[id].begin();
  return am_out_edges[am_out_edges.size()-1].end();
}

adjacent_edges_iterator double_adjacency_map::out_edges_end(id_type id) const
{
  if(is_valid_node(id)) return am_out_edges[id].end();
  return am_out_edges[am_out_edges.size()-1].end();
}

wt_t double_adjacency_map::edge_weight(id_type from_id, id_type to_id) const
{
  if(is_valid_node(from_id) && is_valid_node(to_id))
  {
      adjacent_edges_iterator aeit = am_out_edges[from_id].find(to_id);
      if(aeit != am_out_edges[from_id].end()) return aeit->second;
  }
  return 0;
}              

id_type double_adjacency_map::insert_node()
{
  am_in_edges.push_back(adjacent_edge_sequence());
  am_out_edges.push_back(adjacent_edge_sequence());
  vw_in_degree.push_back(0);
  vw_out_degree.push_back(0);
  return am_in_edges.size();
}

bool double_adjacency_map::insert_edge(id_type from_id, id_type to_id,wt_t weight)
{
  if(!weight) return false;
  if(edge_weight(from_id,to_id)) return false;
  am_in_edges[to_id].insert(make_pair(from_id,weight));
  am_out_edges[from_id].insert(make_pair(to_id,weight));
  st_num_edges++;
  wt_total_wt+= weight;
  vw_in_degree[to_id]+=weight;
  vw_out_degree[from_id]+=weight;
  if(from_id == to_id)
  {
      st_num_self_edges++;
      wt_self_edge_wt+=weight;
  }
  return true;
}

bool double_adjacency_map::delete_edge(id_type from_id,id_type to_id)
{
  wt_t weight = edge_weight(from_id,to_id);
  if(!weight) return false;
  am_in_edges[to_id].erase(from_id);
  am_out_edges[from_id].erase(to_id);
  st_num_edges--;
  wt_total_wt-= weight;
  vw_in_degree[to_id]-=weight;
  vw_out_degree[from_id]-=weight;
  if(from_id == to_id)
  {
      st_num_self_edges--;
      wt_self_edge_wt-=weight;
  }
  return true;
}

wt_t double_adjacency_map::set_edge_wt(id_type from_id,id_type to_id, wt_t weight)
{
  wt_t old_weight = edge_weight(from_id,to_id);
  if(old_weight && weight)
  {
      wt_total_wt += (weight-old_weight);
      vw_in_degree[to_id] += (weight-old_weight);
      vw_out_degree[from_id] += (weight-old_weight);
      if(from_id == to_id) wt_self_edge_wt += (weight-old_weight);
      am_in_edges[to_id][from_id] += (weight-old_weight);
      am_out_edges[from_id][to_id] += (weight-old_weight);
  }
  else 
  {
      if(!old_weight && weight) insert_edge(from_id,to_id,weight);
      else if(old_weight && !weight) delete_edge(from_id,to_id); 
  }
  return old_weight;
}

bool double_adjacency_map::delete_node(id_type id)
{
    if(!is_valid_node(id)) return false;
    id_type last_id = am_in_edges.size()-1;
    deque< pair<id_type,id_type> > edges_to_delete;
    deque< pair< pair<id_type,id_type>,wt_type> >edges_to_add;
    for(adjacent_edges_iterator aeit = am_in_edges[id].begin();aeit != am_in_edges[id].end();aeit++)
        edges_to_delete.push_back(make_pair(aeit->first,id));
    for(adjacent_edges_iterator aeit = am_out_edges[id].begin();aeit != am_out_edges[id].end();aeit++)
        edges_to_delete.push_back(make_pair(id,aeit->first));
    for(adjacent_edges_iterator aeit = am_in_edges[last_id].begin();aeit != am_in_edges[last_id].end();aeit++)
    {
        edges_to_add.push_back(make_pair(make_pair(aeit->first,id),aeit->second));
        edges_to_delete.push_back(make_pair(aeit->first,last_id));
    }
    for(adjacent_edges_iterator aeit = am_out_edges[last_id].begin();aeit != am_out_edges[last_id].end();aeit++)
    {
        edges_to_add.push_back(make_pair(make_pair(id,aeit->first),aeit->second));
        edges_to_delete.push_back(make_pair(last_id,aeit->first));
    }
    for(id_type i=0;i<edges_to_delete.size();i++) delete_edge(edges_to_delete[i].first,edges_to_delete[i].second);
    for(id_type i=0;i<edges_to_add.size();i++) insert_edge(edges_to_add[i].first.first,edges_to_add[i].first.second,edges_to_add[i].second);
    am_in_edges.pop_back();
    am_out_edges.pop_back();
    vw_in_degree.pop_back();
    vw_out_degree.pop_back();
    return true;
}

//
//bool double_adjacency_map::delete_node(id_type id)
//{
//  if(!is_valid_node(id)) return false;
//  for(id_type i=0; i< am_in_edges.size();i++)
//  {
//      if(i != id)
//      {
//        delete_edge(i,id);
//        delete_edge(id,i);
//      }
//      else delete_edge(i,id);
//  }
//  id_type last_id = am_in_edges.size()-1;
//  for(id_type i=0; i< am_in_edges.size();i++)
//  {
//      wt_t weight_in = edge_weight(i,last_id);
//      if(weight_in)
//      {
//          am_out_edges[i].erase(last_id);
//          am_out_edges[i].insert(make_pair(id,weight_in));
//      }
//      wt_t weight_out = edge_weight(last_id,i);
//      if(weight_out)
//      {
//          am_in_edges[i].erase(last_id);
//          am_in_edges[i].insert(make_pair(id,weight_out));
//      }
//  }
//  swap(am_in_edges[id],am_in_edges[last_id]);
//  swap(am_out_edges[id],am_out_edges[last_id]);
//  swap(vw_in_degree[id],vw_in_degree[last_id]);
//  swap(vw_out_degree[id],vw_out_degree[last_id]);
//  am_in_edges.pop_back();
//  am_out_edges.pop_back();
//  vw_in_degree.pop_back();
//  vw_out_degree.pop_back();
//  return true;
//}

bool double_adjacency_map::delete_all_edges()
{
    
    id_type num_nodes = am_in_edges.size();
    if(!num_nodes || !st_num_edges) return false;
    st_num_edges = 0;
    st_num_self_edges = 0;
    wt_total_wt = 0;
    wt_self_edge_wt = 0;
    am_in_edges.assign(num_nodes,adjacent_edge_sequence());
    am_out_edges.assign(num_nodes,adjacent_edge_sequence());
    vw_in_degree.assign(num_nodes,wt_t());
    vw_out_degree.assign(num_nodes,wt_t());
    return true;
}
bool double_adjacency_map::clear()
{
    delete_all_edges();
    id_type num_nodes = am_in_edges.size();
    if(!num_nodes) return false;
    am_in_edges.clear();
    am_out_edges.clear();
    vw_in_degree.clear();
    vw_out_degree.clear();
    return true;
}