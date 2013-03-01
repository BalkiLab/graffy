/* 
 * File:   xmeans_ext.h
 * Author: bharath
 *
 * Created on 1 March, 2013, 10:20 PM
 */

#ifndef XMEANS_EXT_H
#define	XMEANS_EXT_H

#include "includes.h"

namespace XMeans
{
    
    struct xmeans_options{
        const string XMEANS_PATH;
        id_type min_k;
        id_type max_k;
        id_type num_splits;
        id_type max_iters;
        string method;
        bool use_AD;
        id_type max_leaf_size;
        id_type min_box_width;
        double cutoff_factor;
        xmeans_options(id_type min_kv, id_type max_kv, id_type num_splitsv){
            XMEANS_PATH = "/home/bharath/NetBeansProjects/auton/kmeans";
            min_k = min_kv;
            max_kv = max_kv;
            num_splits = num_splitsv;
            use_AD = false;
        }
    };
    
    int run_xmeans_with_file(const string& infile, id_type min_clusters,id_type max_clusters, id_type num_splits, , id_type max_iters,vector<id_type>& labels,vector< vector<double> >& centroids, bool use_ad) {
        const string 
        const string XMEANS_DEFAULT_OPS = ""
        ostringstream oss, oss1, oss2;

        if(use_ad) oss << XMEANS_PATH + "/kmeans kmeans -S AD -k " << min_clusters << " " << "-method blacklist -max_leaf_size 40 -min_box_width 0.03 -cutoff_factor 0.5 -create_universe -max_iter 200 -num_splits 12 -max_ctrs " << max_clusters << " -printclusters " << infile + ".ktemp -printcentroids" << infile + ".ctemp" << " -in " << infile << " >/dev/null 2>/dev/null";
        else oss << XMEANS_PATH + "/kmeans kmeans -k 2 -method blacklist -max_leaf_size 40 -min_box_width 0.03 -cutoff_factor 0.5 -create_universe -max_iter 200 -num_splits 12 -max_ctrs " << max_clusters << " -printclusters " << infile + ".ktemp -printcentroids" << infile + ".ctemp" << " -in " << infile << " >/dev/null 2>/dev/null";
        int r1 = system(oss.str().c_str());
        oss1 << XMEANS_PATH + "/" << "membership " << infile + ".ktemp >" << infile + ".kclusts";
        int r2 = system(oss1.str().c_str());
        labels.clear();
        CDLib::read_vector_1D(infile + ".kclusts",labels);
        CDLib::read_vector_of_vector(infile + ".kclusts",centroids);
        oss2 << "rm " << infile + ".universe" << " " << infile + ".ktemp" << " " << infile + ".kclusts" << " " << infile + ".kc";
        int r3 = system(oss2.str().c_str());
        return r1 && r2 && r3;
    }
    
    void write_xmeans_data_file(const vector< vector<double> >& data, const string& outfile){
        ofstream ofs(outfile);
        for(id_type i=0;i<data[0].size();i++)
            ofs << "x" << i << " ";
        ofs << endl;
        for(id_type i=0;i<data.size();i++)
        {
            for(id_type j=0;j<data[i].size();j++)
                ofs << data[i][j] << " ";
            ofs <<endl;
        }
    }
    
    //data-> columns are features
    void xmeans(const vector< vector<double> >& data, id_type max_k, vector<id_type>& labels, vector<double>& centers,bool use_ad){
        if(!data.empty()){
            string tempfile(tmpnam(NULL));
            write_xmeans_data_file(data,tempfile);
            run_xmeans_with_file(tempfile,maxk,labels,centroids,use_ad);
        }
    }
}


#endif	/* XMEANS_EXT_H */

