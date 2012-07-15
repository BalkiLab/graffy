/* 
 * File:   binary_heap.cpp
 * Author: bharath
 * 
 * Created on 4 April, 2012, 6:32 PM
 */



#include <deque>

#include "binary_heap.h"
using namespace CDLib;

bool binary_heap::compare(const pair<id_type,wt_type>& left,const pair<id_type,wt_type>& right) const { return (b_max) ? (left.second > right.second) : (left.second < right.second) ; }
id_type binary_heap::size() const { return vpiw_heap.size();}
bool binary_heap::empty() const { return !vpiw_heap.size(); }
pair<id_type,wt_type> binary_heap::top() const { return vpiw_heap[0]; }

void binary_heap::heapify_up(id_type pos)
{
     while(pos>0 && compare(vpiw_heap[pos],vpiw_heap[(pos-1)/2]))
     {
            swap(umis_pos[vpiw_heap[pos].first],umis_pos[vpiw_heap[(pos-1)/2].first]);
            //cout << "Heapify_Up("<< pos <<") Swapping (" << vpiw_heap[pos].first << "," << vpiw_heap[pos].second << ") and (" << vpiw_heap[(pos-1)/2].first << "," << vpiw_heap[(pos-1)/2].second << ")" <<endl;
            swap(vpiw_heap[pos],vpiw_heap[(pos-1)/2]);
            pos = (pos-1)/2;
     }
}

void binary_heap::update_key(const pair<id_type,wt_type>& piw_in)
{
    id_type pos = umis_pos[piw_in.first];
    double old_val = vpiw_heap[pos].second;
    vpiw_heap[pos].second = piw_in.second;
    bool logic = b_max ?(old_val < piw_in.second):(old_val > piw_in.second);
    if(logic)heapify_up(pos);
    else heapify_down(pos);
}

void binary_heap::insert(const pair<id_type,wt_type>& piw_in)
{
    vpiw_heap.push_back(piw_in);
    id_type pos = vpiw_heap.size()-1;
    umis_pos.insert(make_pair(piw_in.first,pos));
    update_key(piw_in);
}

void binary_heap::heapify_down(id_type root)
{
    if(root < vpiw_heap.size())
    {
        id_type left = 2*root+1,right = 2*root+2,top = root;
        if(left < vpiw_heap.size() && compare(vpiw_heap[left],vpiw_heap[top])) top = left;
        if(right < vpiw_heap.size() && compare(vpiw_heap[right],vpiw_heap[top])) top = right;
        if(top != root)
        {
            swap(umis_pos[vpiw_heap[top].first],umis_pos[vpiw_heap[root].first]);
            swap(vpiw_heap[top],vpiw_heap[root]);
            heapify_down(top);
        }
    }
}

void binary_heap::pop()
{
    umis_pos.erase(vpiw_heap[0].first);
    umis_pos[vpiw_heap[vpiw_heap.size()-1].first] = 0;    
    vpiw_heap[0] = vpiw_heap[vpiw_heap.size()-1];
    vpiw_heap.pop_back();
    heapify_down(0);
}

binary_heap::binary_heap(const vector<wt_type>& v,bool max): vpiw_heap(),umis_pos(),b_max(max)
{
    for(id_type i=0;i<v.size();i++)
    {
        vpiw_heap.push_back(make_pair(i,v[i]));
        umis_pos.insert(make_pair(i,i));
    }
    for(id_type i=(vpiw_heap.size()+1)/2;i>0;i--) heapify_down(i);
    heapify_down(0);
        
}

binary_heap::binary_heap(bool max) : vpiw_heap(),umis_pos(),b_max(max) {}


