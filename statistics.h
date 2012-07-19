/* 
 * File:   statistics.h
 * Author: bharath
 *
 * Created on May 11, 2012, 10:10 AM
 */

#ifndef STATISTICS_H
#define	STATISTICS_H

#include "typedefs.h"
using namespace std;
namespace CDLib
{

 template <typename T>
    struct statistics
    {
        T min_val;
        T median_val;
        T max_val;
        double mean_val;
        double variance;
    };
    
    template <typename T>
    void compute_statistics(const vector<T>& vec,statistics<T>& res)
    {
        res.variance = 0.0;
        res.mean_val = 0.0;
        res.max_val = (numeric_limits<T>::has_infinity) ? -numeric_limits<T>::infinity() : -numeric_limits<T>::max();
        res.min_val = (numeric_limits<T>::has_infinity) ? numeric_limits<T>::infinity() : numeric_limits<T>::max();
        for(typename vector<T>::const_iterator it = vec.begin();it != vec.end();it++)
        {
            if(res.max_val < *it) res.max_val = *it;
            if(res.min_val > *it) res.min_val = *it;
            res.mean_val += *it;
            res.variance += (*it)*(*it);
        }
        res.mean_val/=(double)vec.size();
        res.variance /= (double)vec.size();
        res.variance -= (res.mean_val*res.mean_val);
        //res.variance /= (double)vec.size();
        vector<T> bleh(vec);
        sort(bleh.begin(),bleh.end());
        res.median_val = (bleh.size()%2) ? (bleh[bleh.size()/2]) :  (bleh[(bleh.size()-1)/2] +  bleh[(bleh.size()+1)/2])/2;
    }
    
    template <typename T>
    string statistics_string(const vector<T>& vec,const string& separator)
    {
        ostringstream oss_ss;
        statistics<T> res;
        compute_statistics<T>(vec,res);
        oss_ss << res.min_val << separator;
        oss_ss << res.median_val << separator; 
        oss_ss << res.max_val << separator;
        oss_ss << res.mean_val << separator;
        oss_ss << res.variance;
        return oss_ss.str();
    }
    
    template <typename T>
    void get_discrete_distribution(const vector<T>& vec,vector<T> & dist)
    {
        statistics<T> stat;
        compute_statistics<T>(vec,stat);
        dist.assign(stat.max_val+1,0);
        for(id_type i=0;i<vec.size();i++) dist[vec[i]]++;
    }

};
#endif	/* STATISTICS_H */

