/* 
 * File:   binary_heap.h
 * Author: bharath
 *
 * Created on 4 April, 2012, 6:32 PM
 */

#ifndef BINARY_HEAP_H
#define	BINARY_HEAP_H
#include "typedefs.h"
using namespace std;
namespace CDLib
{
    class binary_heap {
    private:
        vector< pair<id_type,wt_type> >vpiw_heap;
        unordered_map<id_type,id_type> umis_pos;
        bool b_max;
        bool compare(const pair<id_type,wt_type>& left,const pair<id_type,wt_type>& right) const;
    public:
        
        id_type size() const;
        bool empty() const;
        pair<id_type,wt_type> top() const;
        void update_key(const pair<id_type,wt_type>& piw_in);
        void insert(const pair<id_type,wt_type>& piw_in);
        void heapify_up(id_type pos);
        void heapify_down(id_type root);
        void pop();
        binary_heap(bool min);
        binary_heap(const vector<wt_type>& v,bool max);
    };
}

#endif	/* BINARY_HEAP_H */

