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
        statistics() {
            min_val =0;         median_val = 0;         max_val = 0;
            mean_val = 0;       variance = 0;
        }
    };
    
    template <typename T1, typename T2>
    int kronecker_delta(const T1 sample1, const T2 sample2) { return (sample1 == sample2)?1:0; }
    
    template <typename T>
    void compute_statistics(const vector<T>& vec,statistics<T>& res)
    {
        res.variance = 0.0;
        res.mean_val = 0.0;
        res.max_val = (numeric_limits<T>::has_infinity) ? -numeric_limits<T>::infinity() : -numeric_limits<T>::max();
        res.min_val = (numeric_limits<T>::has_infinity) ? numeric_limits<T>::infinity() : numeric_limits<T>::max();
        res.median_val = 0;
        if(vec.size() > 0){
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
            res.variance /= (double)vec.size();
            vector<T> bleh(vec);
            sort(bleh.begin(),bleh.end());
            res.median_val = (bleh.size()%2) ? (bleh[bleh.size()/2]) :  (bleh[(bleh.size()-1)/2] +  bleh[(bleh.size()+1)/2])/2;
        }
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
    
    template <typename T>
    void get_discrete_distribution_fullvector(const vector<T>& vec,vector<T> & dist)
    {
        dist.assign(vec.size(),0);
        for(id_type i=0;i<vec.size();i++) dist[vec[i]]++;
    }
    
    template <typename T>
    double mean(const vector<T>& samples)
    {
        if (samples.size() == 0)
            return 0;
        if (samples.size() == 1)
            return (double)samples[0];
        T sum = 0;
        sum = accumulate(samples.begin(), samples.end(),0);
        return ((double)sum/samples.size());
    }

    template <typename T>
    double variance(const vector<T>& samples)
    {
        if (samples.size() <= 1)
            return 0;
        double sum_sq=0,sum=0;
        for (id_type i=0 ; i<samples.size(); i++) {
            sum_sq += samples[i] * samples[i];
            sum += samples[i];
        }
        return ((sum_sq/samples.size()) - pow((sum/samples.size()),2));
    }

    template <typename T>
    double std(const vector<T>& samples)
    {
        if (samples.size() <= 1)
            return 0;
        return sqrt(variance(samples));
    }

    template <typename T1, typename T2>
    double covariance(const vector<T1>& sample1, const vector<T2>& sample2)
    {
        if (sample1.size() != sample2.size())
            return -99999;          // Reporting unequal sample size error.
        if (sample1.size() <= 1)
            return 0;
        double mean1 = mean(sample1);
        double mean2 = mean(sample2);
        double sum = 0;
        for (id_type i=0 ; i<sample1.size(); i++)
            sum += (sample1[i] - mean1) * (sample2[i] - mean2);
        return (sum/sample1.size());
    }

    template <typename T>
    void change_to_zero_mean(vector<T>& samples)
    {
        if (samples.size() == 0)
            return;
        if (samples.size() == 1){
            samples[0] = 0;
            return;
        }
        double mean_val = mean(samples);
        for (id_type i=0 ; i<samples.size(); i++)
            samples[i] -= mean_val;
    }
    
    template <typename T>
    void change_to_zero_mean_unit_variance(vector<T>& samples)
    {
        if (samples.size() == 0)
            return;
        if (samples.size() == 1){
            samples[0] = 0;
            return;
        }
        double mean_val = mean(samples);
        double std_val = std(samples);
        for (id_type i=0 ; i<samples.size(); i++)
            samples[i] = (samples[i] - mean_val)/std_val;
    }
    
    template <typename T1, typename T2>
    double pearson_correlation(const vector<T1>& sample1, const vector<T2>& sample2)
    {
        if (sample1.size() != sample2.size())
            return -99999;          // Reporting unequal sample size error.
        if (sample1.size() == 0)
            return 0;
        return covariance(sample1,sample2)/(std(sample1) * std(sample2));
    }
    bool is_distribution(const vector<double>& distribution);
    double distribution_entropy(const vector<double>& distribution);
    double kl_divergence(const vector<double>& distribution1,const vector<double>& distribution2);
    double kl_divergence_symmetric(const vector<double>& distribution1,const vector<double>& distribution2);
    double bhattacharyya_distance(const vector<double>& distribution1,const vector<double>& distribution2);
    double hellinger_distance(const vector<double>& distribution1,const vector<double>& distribution2);

};

#endif	/* STATISTICS_H */

