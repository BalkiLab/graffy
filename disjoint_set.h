/* 
 * File:   disjoint_set.h
 * Author: bharath
 *
 * Created on 12 April, 2012, 10:37 PM
 */

#ifndef DISJOINT_SET_H
#define	DISJOINT_SET_H
#include "typedefs.h"
using namespace std;
namespace CDLib
{
    typedef unordered_map< id_type,pair<id_type,id_type> > id_node_map;
    typedef id_node_map::const_iterator ds_iterator;
    class disjoint_set {
    private:
        id_type st_num_sets;
        id_node_map inm_elems;
    public:
        disjoint_set();
        id_type size() const;
        id_type num_sets() const;
        ds_iterator begin() const;
        ds_iterator end() const;
        pair<ds_iterator,bool> make_set(id_type x);
        pair<ds_iterator,bool> find(id_type x);
        pair<ds_iterator,bool> join(id_type x, id_type y);
    };
}

#endif	/* DISJOINT_SET_H */

