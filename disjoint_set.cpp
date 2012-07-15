/* 
 * File:   disjoint_set.cpp
 * Author: bharath
 * 
 * Created on 12 April, 2012, 10:37 PM
 */

#include "disjoint_set.h"
using namespace CDLib;

disjoint_set::disjoint_set() : st_num_sets(0),  inm_elems() {}

id_type disjoint_set::size() const { return inm_elems.size(); }
id_type disjoint_set::num_sets() const { return st_num_sets; }

ds_iterator disjoint_set::begin() const { return inm_elems.begin(); }
ds_iterator disjoint_set::end() const { return inm_elems.end(); }

pair<ds_iterator,bool> disjoint_set::make_set(id_type x) { return inm_elems.insert(make_pair(x,make_pair(x,0))); }

pair<ds_iterator,bool> disjoint_set::find(id_type x)
{
    ds_iterator dit = inm_elems.find(x) ;
    if(dit ==  inm_elems.end()) return make_pair(dit,false);
    if(inm_elems[x].first != x)
        inm_elems[x].first = find(inm_elems[x].first).first->second.first;
    return make_pair(inm_elems.find(dit->second.first),true);
}

pair<ds_iterator,bool> disjoint_set::join(id_type x, id_type y)
{
    pair<ds_iterator,bool> ret_x = find(x);
    pair<ds_iterator,bool> ret_y = find(y);
    if(!ret_x.second || !ret_y.second ) return make_pair(inm_elems.end(),false);
    if(ret_x.first == ret_y.first) return ret_x;
    if(ret_x.first->second.second < ret_y.first->second.second)
    {
        inm_elems[ret_x.first->first].first = ret_y.first->first;
        return ret_y;
    }
    else
    {
        inm_elems[ret_y.first->first].second++;
        inm_elems[ret_y.first->first].first = ret_x.first->first;
        return ret_x;
    }
}
