/* 
* File:   bidirectional_label_map.h
* Author: bharath
*
* Created on April 1, 2012, 11:00 AM
*/

#ifndef BIDIRECTIONAL_LABEL_MAP_H
#define	BIDIRECTIONAL_LABEL_MAP_H
#include <unordered_map>
#include <string>
#include <sstream>
#include "typedefs.h"
using namespace std;
namespace CDLib
{ 
  typedef unordered_map<string,id_type> forward_map;
  typedef unordered_map<id_type,string> reverse_map;
  typedef forward_map::const_iterator node_label_iterator;
  class bidirectional_label_map 
  {
  private:
      forward_map fm_labels;
      reverse_map rm_ids;
  public:
      
       id_type size() const;
       node_label_iterator begin() const;
       node_label_iterator end() const;
       string get_label(id_type id) const;
       id_type get_id(const string& label) const;
       bool insert(const string& label);
       bool erase(const string& label);
       bool erase(id_type id);
       bool clear();
  };
};

#endif	/* BIDIRECTIONAL_LABEL_MAP_H */

