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

#endif	/* UTITLITY_H */

