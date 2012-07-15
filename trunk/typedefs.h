/* 
 * File:   typedefs.h
 * Author: bharath
 *
 * Created on 14 April, 2012, 6:14 PM
 */

#ifndef TYPEDEFS_H
#define	TYPEDEFS_H

#include <algorithm>
#include <deque>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <stack>
#include <map>
#include <iterator>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <numeric>
#include <limits>
#include <list>
#include <set>
#include <sstream>
#include <string>
#include <queue>
#include <vector>
#include <omp.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <utility>

namespace std {
    template <class T1, class T2>
        class hash< pair<T1,T2> >{
        public :
            size_t operator()(const pair<T1,T2>& my_pair ) const
            {
                return hash<T1>()(my_pair.first) ^ hash<T2>()(my_pair.second);
            }
    };
};

#define BILLION  1000000000L;
using namespace std;
namespace CDLib
{
    //Basic Types
    typedef unsigned long id_type;
    typedef double wt_type;
    typedef double wt_t;

    //Paths and Components, Community Detection Algorithms
    typedef unordered_set<id_type> node_set;
    typedef unordered_set<string> node_set_string;
    
};






#endif	/* TYPEDEFS_H */

