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
        string XMEANS_PATH;
        id_type min_k;
        id_type max_k;
        id_type num_splits;
        id_type max_iters;
        string method;
        bool use_AD;
        double max_leaf_size;
        double min_box_width;
        double cutoff_factor;
        void default_options(){
            XMEANS_PATH = "/opt/Extensions/auton/kmeans";
            min_k = 2;
            max_k = 2;
            num_splits = 24;
            use_AD = false;
            method = "blacklist";
            max_leaf_size = 40;
            min_box_width = 0.03;
            cutoff_factor = 0.5;
            max_iters = 100;
        }
        xmeans_options(id_type min_kv, id_type max_kv){
            default_options();
            min_k = min_kv;
            max_k = max_kv;
        }
        xmeans_options(){
            default_options();
        }
        int run_xmeans_command(const string& fileprefix,vector<id_type>& labels,vector< vector<double> >& centroids) const{
            ostringstream oss,ossm,ossr;
            oss << XMEANS_PATH << "/kmeans kmeans ";
            if(use_AD) oss << "-S AD ";
            oss << "-k " << min_k << " ";
            oss << "-method " << method << " ";
            oss << "-max_leaf_size " << max_leaf_size << " ";
            oss << "-min_box_width " << min_box_width << " ";
            oss << "-cutoff_factor " << cutoff_factor << " ";
            oss << "-create_universe ";
            oss << "-num_splits " << num_splits << " ";
            oss << "-max_ctrs " << max_k << " ";
            oss << "-max_iters " << max_iters << " ";
            oss << "-printclusters " << fileprefix + ".ktemp ";
            oss << "-printcentroids " << fileprefix + ".ctemp ";
            oss << "-in " << fileprefix << " ";
            oss << "2>"<< fileprefix + ".log >" << fileprefix + ".out";
            int r1 = system(oss.str().c_str());
            ossm << XMEANS_PATH << "/membership " << fileprefix << ".ktemp >" << fileprefix << ".kclusts";
            int r2 = system(ossm.str().c_str());
            read_vector_1D(fileprefix + ".kclusts",labels);
            read_vector_of_vector(fileprefix + ".ctemp",centroids);
            ossr << "rm " << fileprefix + ".universe" << " " << fileprefix + ".ktemp" << " " << fileprefix + ".kclusts" << " " << fileprefix + ".ctemp" << " " << fileprefix + ".log"  << " " << fileprefix + ".out";
            int r3 = system(ossr.str().c_str());
            return r1 && r2 && r3;
        }

    };


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
    void compute_centroids(const vector< vector<double> >& data, vector<id_type>& labels, vector <vector<double> >& centroids){
        id_type num_comms = reindex_communities(labels);
        centroids.assign(num_comms,vector<double>(data[0].size(),0));
        vector<id_type> counts(num_comms,0);
        for(id_type i=0;i<data.size();i++){
            counts[labels[i]]++;
            for(id_type j=0;j<data[0].size();j++)
                centroids[labels[i]][j] += data[i][j];
        }
        for(id_type i=0;i<num_comms;i++)
            for(id_type j=0;j<data[0].size();j++)
                centroids[i][j]/= static_cast<double>(counts[i]);
    }

    //data-> columns are features
    void xmeans(const vector< vector<double> >& data, const xmeans_options& opts,vector<id_type>& labels, vector<vector<double> >& centroids){
        if(!data.empty()){
            string tempfile(tmpnam(NULL));
            write_xmeans_data_file(data,tempfile);
            opts.run_xmeans_command(tempfile,labels,centroids);
            compute_centroids(data,labels,centroids);
        }
    }
}


//xmeans(data,xmeans_options(),labels,centroids);
//xmeans(data,xmeans_options(min_k,max_k),labels,centroids);
//xopt.sds = adad
//xmeans(data,xopt,labels,centroids);

#endif	/* XMEANS_EXT_H */

