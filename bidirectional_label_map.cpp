/* 
* File:   bidirectional_label_map.cpp
* Author: bharath
* 
* Created on April 1, 2012, 11:00 AM
*/

#include "bidirectional_label_map.h"
using namespace CDLib;

id_type bidirectional_label_map::size() const { return fm_labels.size(); } 
      
node_label_iterator bidirectional_label_map::begin() const { return fm_labels.begin(); } 

node_label_iterator bidirectional_label_map::end() const { return fm_labels.end(); } 

string bidirectional_label_map::get_label(id_type id) const
{
  unordered_map<id_type,string>::const_iterator it = rm_ids.find(id);
  if(it != rm_ids.end()) return it->second;
  return string("");
}

id_type bidirectional_label_map::get_id(const string& label) const
{
  node_label_iterator it = fm_labels.find(label);
  if(it != end()) return it->second;
  return size();
}

bool bidirectional_label_map::insert(const string& label)
{
  id_type curr_size = size();
  if(get_id(label) == curr_size)
  {
      fm_labels.insert(make_pair(label,curr_size));
      rm_ids.insert(make_pair(curr_size,label));
      return true;
  }
  return false;    
}

bool bidirectional_label_map::erase(const string& label)
{
  id_type del_id = get_id(label);
  if(del_id != size())
  {
      fm_labels.erase(label);
      rm_ids.erase(del_id);
      return true;
  }
  return false;
}

bool bidirectional_label_map::erase(id_type id)
{
  id_type last_id = size() -1;
  string del_label = get_label(id),rep_label = get_label(last_id);
  if(del_label != "")
  {
      fm_labels[rep_label] = id; 
      rm_ids[id] = rep_label;
      fm_labels.erase(del_label);
      rm_ids.erase(last_id);
      return true;
  }
  return false;
}

bool bidirectional_label_map::clear()
{
    if(!fm_labels.size()) return false;
    fm_labels.clear();
    rm_ids.clear();
    return true;
}
