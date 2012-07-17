/* 
 * File:   utitlity.h
 * Author: bharath
 *
 * Created on May 13, 2012, 6:04 AM
 */

#ifndef UTITLITY_H
#define	UTITLITY_H

#include "typedefs.h"
#include "random.h"
using namespace CDLib;

template <typename T>
T str2T(const string& str)
{
    T retval;
    istringstream iss(str);
    iss >> retval;
    return retval;
}

template <typename T>
string T2str(T val)
{
    ostringstream oss;
    oss << val;
    return oss.str();
}

template <typename T>
vector<T> split(const string& s,const char delim)
{
    vector<T> elems;
    stringstream ss(s);
    string item;
    while(getline(ss,item,delim))
        elems.push_back(str2T<T>(item));
    return elems;
}

template <typename T>
void cumsum(const vector<T>& v_in,vector<T>& v_out)
{
    v_out.assign(v_in.size(),0);
    v_out[0] = v_in[0];
    for(id_type i=1;i<v_in.size();i++)
        v_out[i] = v_out[i-1] + v_in[i];
}

template <class Element, class Val>
class max_manager
{
private:
    deque<Element> max_container;
    Val max_val;
    unsigned long cont_last_index;
public:
    typedef typename deque<Element>::const_iterator max_iterator;
    max_manager() : max_container(), max_val(),cont_last_index(0)
    {
        max_val = (numeric_limits<Val>::has_infinity) ? -numeric_limits<Val>::infinity() : -numeric_limits<Val>::max();
    }
    id_type size() const { return cont_last_index+1; }
    Element get_max() const 
    {
        RandomGenerator<size_t> gen(0,cont_last_index,1);
        return max_container[gen.next()]; 
    }
    Val get_max_val() const { return max_val; }
    max_iterator begin() { return max_container.begin(); }
    max_iterator end() { return begin()+ cont_last_index+1;}
    bool insert(Element elem, Val m_val)
    {
        if(max_val > m_val) return false;
        else if(max_val == m_val)
        {
            max_container.push_back(elem);
            cont_last_index = ((max_container.empty()) ? 0 : cont_last_index+1);
        }
        else
        {
            if(max_container.empty()) max_container.push_back(elem);
            else max_container[0] = elem;
            cont_last_index =0;
        }
        return true;
    }
};

    template <class Element, typename Priority>
    struct min_heap
    {
       bool operator()(const pair<Element,Priority>& left,const pair<Element,Priority>& right) const
       {
          return (left.second < right.second) ;
       }
    };
    template <class Element, typename Priority>
    struct max_heap
    {
       bool operator()(const pair<Element,Priority>& left,const pair<Element,Priority>& right) const
       {
          return (left.second > right.second) ;
       }
    };
    
    
    
    template <class Element, class Priority,class Compare>
    class hash_heap
    {
    private:
        deque< pair<Element,Priority> >vpiw_heap;
        unordered_map<Element,id_type> umis_pos;
        Compare cmp;
    public:
        typedef typename deque< pair<Element,Priority> >::const_iterator iterator;
        id_type size() const { return vpiw_heap.size();}
        bool empty() const { return !vpiw_heap.size(); }
        pair<Element,Priority> front() const { return vpiw_heap[0]; }
        iterator begin() const { return vpiw_heap.begin(); }
        iterator end() const { return vpiw_heap.end(); }
        iterator find(Element elem) const 
        { 
            typename unordered_map<Element,id_type>::const_iterator it = umis_pos.find(elem);
            if(it != umis_pos.end()) return vpiw_heap.begin() + it->second;
            else return end();
        }
        bool update_key(const pair<Element,Priority>& piw_in)
        {
            if(umis_pos.find(piw_in.first) != umis_pos.end())
            {
                id_type pos = umis_pos[piw_in.first];
                pair<Element,Priority> old_val = vpiw_heap[pos];
                vpiw_heap[pos].second = piw_in.second;
                if(!cmp(piw_in,old_val))heapify_up(pos);
                else heapify_down(pos);
                return true;
            }
            return insert(piw_in);
        }
        bool insert(const pair<Element,wt_type>& piw_in)
        {
            if(umis_pos.find(piw_in.first) != umis_pos.end()) return update_key(piw_in);
            vpiw_heap.push_back(piw_in);
            id_type pos = vpiw_heap.size()-1;
            umis_pos.insert(make_pair(piw_in.first,pos));
            update_key(piw_in);
            return true;
        }
        void heapify_up(id_type pos)
        {
            while(pos>0 && cmp(vpiw_heap[pos],vpiw_heap[(pos-1)/2]))
            {
                    swap(umis_pos[vpiw_heap[pos].first],umis_pos[vpiw_heap[(pos-1)/2].first]);
                    swap(vpiw_heap[pos],vpiw_heap[(pos-1)/2]);
                    pos = (pos-1)/2;
            }
        }
        void heapify_down(id_type root)
        {
            if(root < vpiw_heap.size())
            {
                id_type left = 2*root+1,right = 2*root+2,top = root;
                if(left < vpiw_heap.size() && cmp(vpiw_heap[left],vpiw_heap[top])) top = left;
                if(right < vpiw_heap.size() && cmp(vpiw_heap[right],vpiw_heap[top])) top = right;
                if(top != root)
                {
                    swap(umis_pos[vpiw_heap[top].first],umis_pos[vpiw_heap[root].first]);
                    swap(vpiw_heap[top],vpiw_heap[root]);
                    heapify_down(top);
                }
            }
        }
        void pop()
        {
            umis_pos.erase(vpiw_heap[0].first);
            umis_pos[vpiw_heap[vpiw_heap.size()-1].first] = 0;    
            vpiw_heap[0] = vpiw_heap[vpiw_heap.size()-1];
            vpiw_heap.pop_back();
            heapify_down(0);
        }
        bool remove(Element elem)
        {
            if(umis_pos.find(elem) == umis_pos.end()) return false;
            id_type pos = umis_pos[elem];
            umis_pos.erase(vpiw_heap[pos].first);
            umis_pos[vpiw_heap[vpiw_heap.size()-1].first] = pos;    
            vpiw_heap[pos] = vpiw_heap[vpiw_heap.size()-1];
            vpiw_heap.pop_back();
            heapify_down(pos);
        }
        hash_heap():vpiw_heap(), umis_pos(),cmp() {}
        hash_heap(const vector< pair<Element,Priority> >& v)
        {
            for(id_type i=0;i<v.size();i++)
            {
                vpiw_heap.push_back(make_pair(v[i].first,v[i].second));
                umis_pos.insert(make_pair(v[i].first,i));
            }
            for(id_type i=(vpiw_heap.size()+1)/2;i>0;i--) heapify_down(i);
            heapify_down(0);
        }
    };
    
    typedef hash_heap<id_type,double,min_heap<id_type,double> > id_dbl_min_heap;
    typedef hash_heap<id_type,double,max_heap<id_type,double> > id_dbl_max_heap;
    

#endif	/* UTITLITY_H */

