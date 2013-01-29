/* 
 * File:   community_tools.cpp
 * Author: bharath
 * 
 * Created on 5 April, 2012, 2:27 AM
 */

#include "community_tools.h"

using namespace CDLib;

CDLib::cluster_edges::cluster_edges() : 
    num_inter_cluster_edges(0),
    num_intra_cluster_edges(0),
    num_expected_inter_cluster_edges(0.0),
    num_expected_intra_cluster_edges(0.0),
    wt_inter_cluster_edges(0.0),
    wt_intra_cluster_edges(0.0),
    wt_expected_inter_cluster_edges(0.0),
    wt_expected_intra_cluster_edges(0.0),
    max_odf(-numeric_limits<double>::infinity()),
    avg_odf(0.0),
    flake_odf(0.0),
    is_weak_radicchi_community(false),
    is_strong_radicchi_community(false)
{}

CDLib::cluster_edges::cluster_edges(const graph& g,node_set& comm)
{
    num_inter_cluster_edges = 0;
    num_intra_cluster_edges = 0;
    num_expected_inter_cluster_edges = 0;
    num_expected_intra_cluster_edges = 0;
    wt_inter_cluster_edges = 0.0;
    wt_intra_cluster_edges = 0.0;
    wt_expected_inter_cluster_edges = 0.0;
    wt_expected_intra_cluster_edges = 0.0;
    is_strong_radicchi_community = true;
    max_odf = -numeric_limits<double>::infinity();
    avg_odf = 0.0;
    flake_odf = 0.0;
    for(node_set::iterator it = comm.begin(); it != comm.end(); it++)
    {
        double node_intra_cluster_edges = 0.0;
        double node_inter_cluster_edges = 0.0;
        for(adjacent_edges_iterator aeit = g.out_edges_begin(*it); aeit != g.out_edges_end(*it); aeit++)
        {
            double wt_expected = (g.get_node_out_weight(*it)*g.get_node_in_weight(aeit->first))/(2*g.get_total_weight());
            double num_expected = (double)(g.get_node_out_degree(*it)*g.get_node_in_degree(aeit->first)/(double)(2*g.get_num_edges()));
            if(comm.find(aeit->first) != comm.end() ) 
            {
                node_intra_cluster_edges+=aeit->second;
                num_intra_cluster_edges++;
                num_expected_intra_cluster_edges += num_expected;
                wt_intra_cluster_edges += aeit->second;
                wt_expected_intra_cluster_edges += wt_expected;
            }
            else
            {
                node_inter_cluster_edges += aeit->second;
                num_inter_cluster_edges++;
                num_expected_inter_cluster_edges += num_expected;
                wt_inter_cluster_edges += aeit->second;
                wt_expected_inter_cluster_edges += wt_expected;
            }
            
        }
        if(node_inter_cluster_edges > node_intra_cluster_edges) is_strong_radicchi_community = false;
        double odf = node_inter_cluster_edges/g.get_node_out_degree(*it);
        max_odf = (odf > max_odf) ? odf : max_odf;
        avg_odf += odf;
        if(node_intra_cluster_edges < (g.get_node_out_degree(*it)/2)) flake_odf++;
    }
    is_weak_radicchi_community = (wt_intra_cluster_edges > wt_inter_cluster_edges);
    avg_odf /= comm.size();
    flake_odf /= (double)comm.size();
}

double CDLib::volume_comm(const graph& g, node_set& comm)
{
    cluster_edges ce(g,comm);
    return ce.wt_intra_cluster_edges;
}

double CDLib::cut_comm(const graph& g, node_set& comm)
{
    cluster_edges ce(g,comm);
    return ce.wt_inter_cluster_edges;
}

double CDLib::ratio_cut_comm(const graph& g, node_set& comm)
{
    cluster_edges ce(g,comm);
    return (double)ce.wt_inter_cluster_edges / (double)(comm.size()*(g.get_num_nodes() - comm.size()));
}

double CDLib::ratio_assoc_comm(const graph& g, node_set& comm)
{
    cluster_edges ce(g,comm);
    return (double)ce.wt_intra_cluster_edges / (double)(comm.size());
}

double CDLib::normalized_assoc_comm(const graph& g, node_set& comm)
{
    cluster_edges ce(g,comm);
    return ce.wt_intra_cluster_edges/(2*ce.wt_intra_cluster_edges + ce.wt_inter_cluster_edges);
}

double CDLib::conductance_comm(const graph& g, node_set& comm)
{
    cluster_edges ce(g,comm);
    return ce.wt_inter_cluster_edges/(2*ce.wt_intra_cluster_edges + ce.wt_inter_cluster_edges);
}

double CDLib::resistance_comm(const graph& g, node_set& comm)
{
    cluster_edges ce(g,comm);
    return ce.wt_inter_cluster_edges/(2*(g.get_total_weight() - ce.wt_intra_cluster_edges) + ce.wt_inter_cluster_edges);
}


double CDLib::expansion_comm(const graph& g, node_set& comm)
{
    cluster_edges ce(g,comm);
    return (double)ce.num_inter_cluster_edges/(double)comm.size();
}

double CDLib::normalized_cut_comm(const graph& g, node_set& comm)
{
    return conductance_comm(g,comm) + resistance_comm(g,comm);
}

double CDLib::internal_density_comm(const graph& g, node_set& comm)
{
    cluster_edges ce(g,comm);
    return ((g.is_directed()? 1 : 2))*ce.num_intra_cluster_edges/(double)(comm.size()*(comm.size()+((g.get_num_self_edges())?-1:1)));
}

double CDLib::max_odf_comm(const graph& g, node_set& comm)
{
    cluster_edges ce(g,comm);
    return ce.max_odf;
}

double CDLib::avg_odf_comm(const graph& g, node_set& comm)
{
    cluster_edges ce(g,comm);
    return ce.avg_odf;
}

double CDLib::flake_odf_comm(const graph& g, node_set& comm)
{
    cluster_edges ce(g,comm);
    return ce.flake_odf;
}

bool CDLib::radicchi_community(const graph& g, node_set& comm,bool strong)
{
    cluster_edges ce(g,comm);
    if(strong) return ce.is_strong_radicchi_community;
    return ce.is_weak_radicchi_community;
}

double CDLib::modularity_comm(const graph& g, node_set& comm)
{
    cluster_edges ce(g,comm);
    return (1/(2*g.get_total_weight()))*(ce.wt_intra_cluster_edges - ce.wt_expected_intra_cluster_edges);
}

double CDLib::modularity_density_comm(const graph& g, node_set& comm)
{
     cluster_edges ce(g,comm);
     return ce.wt_intra_cluster_edges/ce.wt_expected_intra_cluster_edges;
}

double CDLib::partition_quality(const graph& g,vector<node_set>& comms,double (*func)(const graph& g,node_set& comm))
{
    double mod_val =0;
    for(id_type i=0;i<comms.size();i++)
        mod_val += func(g,comms[i]);
    return mod_val;  
}

double CDLib::modularity(const graph& g, vector<node_set>& comms)
{
    double mod_val =0;
    for(id_type i=0;i<comms.size();i++)
    {
        cluster_edges ce(g,comms[i]);
        mod_val += (ce.wt_intra_cluster_edges - ce.wt_expected_intra_cluster_edges);
    }
    return mod_val/(2*g.get_total_weight());      
}

double CDLib::modularity_density(const graph& g, vector<node_set>& comms)
{
    double mod_val;
    for(id_type i=0;i<comms.size();i++)
    {
        cluster_edges ce(g,comms[i]);
        mod_val += (ce.wt_intra_cluster_edges/ce.wt_expected_intra_cluster_edges);
    }
    return mod_val;      
}

double CDLib::community_score(const graph& g, vector<node_set>& comms)
{
    double mod_val;
    for(id_type i=0;i<comms.size();i++)
    {
        cluster_edges ce(g,comms[i]);
        mod_val += pow((2*ce.wt_intra_cluster_edges)/comms[i].size(),2);
    }
    return mod_val;      
}

double npr(id_type n, id_type r)
{
    double ret_val = 1;
    for(id_type i=0;i<r;i++)ret_val *= (n-i);
    return ret_val;
}
double ncr(id_type n, id_type r) 
{ 
    double ret_val_num=npr(n,r),ret_val_denom = 1;
    for(id_type i=0;i<r;i++)ret_val_denom *= i;
    return ret_val_num/ret_val_denom;
}

double CDLib::description_length(const graph& g, vector<node_set>& comms)
{
    if(!comms.size() || !g.get_num_edges()) return 0;
    double mod_val = g.get_num_nodes()*log(comms.size()) + 0.5*comms.size()*(comms.size()+1)*log(g.get_num_edges());
    id_type prod1 = 1,prod2 = 1;
    for(id_type i=0; i<comms.size();i++)
    {
        cluster_edges ce(g,comms[i]);
        prod1 *= ncr((comms[i].size()*(comms[i].size()-1)/2),ce.num_intra_cluster_edges);
        for(id_type j=0;j<i;j++)
        {
            id_type cij = 0;
            for(node_set::iterator it = comms[i].begin(); it != comms[i].end(); it++)
                for(adjacent_edges_iterator aeit = g.out_edges_begin(*it);aeit != g.out_edges_end(*it);aeit++)
                    if(comms[j].find(aeit->first)!= comms[j].end())
                        cij++;
            prod2 *= ncr(comms[i].size()*comms[j].size(),cij);
        }
    }
    cerr <<"DEBUG "<< prod1 << " " << prod2 << " "  << mod_val<< endl;
    return mod_val + (prod1 && prod2) ? log(prod1*prod2) : 0;
}


void CDLib::compute_imp_metrics_partition(const graph& g, vector<node_set>& comms,vector<double>& metrics)
{
    metrics.assign(4,0.0);
    for(id_type i=0;i<comms.size();i++)
    {
        cluster_edges ce(g,comms[i]);
        metrics[0] +=  ce.wt_intra_cluster_edges/static_cast<double>(comms[i].size()); // RASSOC
        metrics[1] +=  ce.wt_inter_cluster_edges/static_cast<double>(comms[i].size()); // EXP
        metrics[2] +=  ce.wt_intra_cluster_edges/(2*ce.wt_intra_cluster_edges + ce.wt_inter_cluster_edges); //NASSOC
        metrics[3] +=  ce.wt_inter_cluster_edges/(2*ce.wt_intra_cluster_edges + ce.wt_inter_cluster_edges); //Conductance
    }
}

void CDLib::compute_all_metrics_partition(const graph& g, vector<node_set>& comms,vector<double>& metrics)
{
    metrics.assign(14,0.0);
    for(id_type i=0;i<comms.size();i++)
    {
        cluster_edges ce(g,comms[i]);
        metrics[0] +=  ce.wt_intra_cluster_edges; // Community Volume
        metrics[1] += ce.wt_inter_cluster_edges; //Cut
        metrics[2] += (double)ce.wt_inter_cluster_edges / (double)(comms[i].size()*(g.get_num_nodes() - comms[i].size())); //Cut Ratio
        metrics[3] +=  ce.wt_inter_cluster_edges/(2*ce.wt_intra_cluster_edges + ce.wt_inter_cluster_edges); //Conductance
        metrics[4] += ce.wt_inter_cluster_edges/(2*(g.get_total_weight() - ce.wt_intra_cluster_edges) + ce.wt_inter_cluster_edges); //Resistance
        metrics[5] += (double)ce.num_inter_cluster_edges/(double)comms[i].size(); //Expansion
        metrics[6] += ((g.is_directed()? 1 : 2))*ce.num_inter_cluster_edges/(double)(comms[i].size()*(comms[i].size()+((g.get_num_self_edges())?-1:1))); //Internal Density
        metrics[7] +=  (ce.wt_inter_cluster_edges/(2*ce.wt_intra_cluster_edges + ce.wt_inter_cluster_edges)) + (ce.wt_inter_cluster_edges/(2*(g.get_total_weight() - ce.wt_intra_cluster_edges) + ce.wt_inter_cluster_edges));//Normalized Cut
        metrics[8] += ce.max_odf/comms[i].size();
        metrics[9] += ce.avg_odf/comms[i].size();
        metrics[10] += ce.flake_odf/comms[i].size();
        metrics[11] += (1/4*(g.get_total_weight())*(ce.wt_intra_cluster_edges - ce.wt_expected_intra_cluster_edges)); //Modularity
        metrics[12] += (ce.wt_intra_cluster_edges/ce.wt_expected_intra_cluster_edges); //Modularity Density
        metrics[13] += pow((2*ce.wt_intra_cluster_edges)/comms[i].size(),2);        
    }
    metrics[11]/=(4*g.get_total_weight());  
}

void CDLib::compute_community_metrics(const graph& g, node_set& comm,community_metrics& metrics)
{
    metrics.size = comm.size();
    metrics.intracluster_edges = 0;
    metrics.intercluster_edges = 0;
    metrics.modularity = 0;
    metrics.modularity_density = 0;
    for(node_set::iterator nit = comm.begin();nit != comm.end();nit++)
    {
        for(adjacent_edges_iterator aeit = g.out_edges_begin(*nit);aeit != g.out_edges_end(*nit);aeit++)
        {
            if(is_member_of(aeit->first,comm))
            {
                metrics.intracluster_edges+=aeit->second;
                metrics.modularity += aeit->second - ((g.get_node_out_weight(*nit)*g.get_node_in_weight(aeit->first))/2*g.get_total_weight());
                metrics.modularity_density += (2*g.get_total_weight()*aeit->second)/(g.get_node_out_weight(*nit)*g.get_node_in_weight(aeit->first));
            }
            else metrics.intercluster_edges+=aeit->second;
        }
    }
    metrics.intracluster_edges/=2;
    double volume = 2*metrics.intracluster_edges + metrics.intercluster_edges;
    metrics.expansion = metrics.intercluster_edges/static_cast<double>(metrics.size);
    metrics.conductance = metrics.intercluster_edges/volume;
    metrics.modularity /= g.get_total_weight();
    metrics.nassoc = metrics.intracluster_edges/volume;
    metrics.cohesion = metrics.intracluster_edges/static_cast<double>(metrics.size);
    metrics.internal_density = ((g.is_directed() ? 2 : 1)* metrics.intracluster_edges)/static_cast<double>(metrics.size*metrics.size-1);
    metrics.community_score = pow(metrics.intracluster_edges,2)/static_cast<double>(metrics.size);
    graph cg(0,0);
    extract_subgraph(g,comm,cg);
    vector<id_type> out_degrees;
    statistics<id_type> out_degree_stats;
    out_degrees.reserve(cg.get_num_nodes());
    for(id_type i=0;i<cg.get_num_nodes();i++) out_degrees.push_back(g.get_node_out_degree(i))   ;
    compute_statistics(out_degrees,out_degree_stats);
    metrics.degree_homogenity = out_degree_stats.variance/out_degree_stats.mean_val;
    double denom = out_degree_stats.mean_val*comm.size();
    metrics.degree_entropy = 0;
    for(id_type i=0;i<out_degrees.size();i++) if(out_degrees[i])metrics.degree_entropy +=(out_degrees[i]/denom)*log(out_degrees[i]/denom);
    vector<double> qinit(cg.get_num_nodes(),1/static_cast<double>(cg.get_num_nodes())),qfinal;
    for(id_type i=0;i<static_cast<id_type>(log(cg.get_num_nodes())/log(2));i++)
    {
        multiply_vector_transform(cg,qinit,transform_func_row_stochastic,qfinal);
        qinit = qfinal;
    }
    qinit.assign(cg.get_num_nodes(),1/static_cast<double>(cg.get_num_nodes()));
    metrics.rwalk_entropy = 0;
    if(cg.get_num_nodes()>=2)for(id_type i=0;i<cg.get_num_nodes();i++) if(qfinal[i] && qinit[i]) metrics.rwalk_entropy += qinit[i]*log(qinit[i]/qfinal[i]);
}

bool comm_find_pair(vector<node_set>& comms,id_type x,id_type y)
{
    for(id_type i=0; i< comms.size();i++)
        if(comms[i].find(x)!= comms[i].end() && comms[i].find(y)!= comms[i].end())
            return true;
    return false;
}

id_type num_common(node_set& comm1, node_set& comm2)
{
    id_type common = 0;
    for(node_set::iterator nsit = comm1.begin(); nsit != comm1.end(); nsit++)
        if(comm2.find(*nsit) != comm2.end()) common++;
    return common;
}

id_type max_common(vector<node_set>& comms1, vector<node_set>& comms2)
{
    id_type sum_max = 0;
    for(id_type i=0; i< comms1.size(); i++)
    {
        id_type max_common = 0;
        for(id_type j=0; j<comms2.size();j++)
        {
            id_type n_common = num_common(comms1[i],comms2[j]);
            max_common = (n_common > max_common) ? n_common : max_common;
        }
        sum_max+= max_common;
    }
    return sum_max;
}

double CDLib::rand_index(id_type num_nodes, vector<node_set>& comms1, vector<node_set>& comms2)
{
    id_type a00=0,a11=0,a01=0,a10=0;
    for(id_type i=0; i<num_nodes; i++)
    {
        for(id_type j=0; j<num_nodes; j++)
        {
            bool find_1 = comm_find_pair(comms1,i,j);
            bool find_2 = comm_find_pair(comms2,i,j);
            if(find_1 && find_2) a11++;
            else if(!find_1 && !find_2) a00++;
            else if(!find_1 && find_2) a01++;
            else a10++;
        }
    }
    return (double)(a00+a11)/(double)(a10 + a01 + a00 + a11);
}

double CDLib::dongen_index(id_type num_nodes, vector<node_set>& comms1, vector<node_set>& comms2)
{
    id_type t1 = max_common(comms1,comms2),t2 = max_common(comms2,comms1);
    return 1 - ((double)(t1+t2)/(double)(2*num_nodes));
}

double CDLib::nmi(id_type num_nodes, vector<node_set>& comms1, vector<node_set>& comms2)
{
    if(!num_nodes || !comms1.size() || !comms2.size()) return 0;
    double num=0,denom1=0,denom2=0;
    for(id_type i=0; i< comms1.size(); i++)    
        if(comms1[i].size())denom1 += ((double)comms1[i].size()/(double)num_nodes)*log((double)comms1[i].size()/(double)num_nodes);
    for(id_type i=0; i< comms2.size(); i++)
        if(comms2[i].size())denom2 += ((double)comms2[i].size()/(double)num_nodes)*log((double)comms2[i].size()/(double)num_nodes);
    for(id_type i=0; i< comms1.size(); i++)
    {
        for(id_type j=0; j<comms2.size();j++)
        {
            id_type n_com = num_common(comms1[i],comms2[j]);
            if(n_com && comms1[i].size() && comms2[j].size()) num += ((double)n_com/(double)num_nodes)*log((double)(n_com*num_nodes)/(double)(comms1[i].size()*comms2[j].size()));
        }
    }
    return (double)(-(2*num))/(double)(denom1+denom2);
}

double CDLib::variation_of_information(id_type num_nodes,vector<node_set>& comms1, vector<node_set>& comms2)
{
    double num1=0,num2=0;
    for(id_type i=0; i< comms1.size(); i++)
    {
        for(id_type j=0; j<comms2.size();j++)
        {
            id_type n_com = num_common(comms1[i],comms2[j]);
            if(comms1[i].size() && n_com && num_nodes) num1 += ((double)n_com/(double)num_nodes)*log((double)n_com/(double)comms1[i].size());
            if(comms2[j].size() && n_com && num_nodes) num2 += ((double)n_com/(double)num_nodes)*log((double)n_com/(double)comms2[j].size());
        }
    }
    return -(num1 + num2);
}

void CDLib::convert_labels_to_communities(vector<id_type>& labels,vector<node_set>& communities)
{
    communities.clear();
    unordered_map<id_type,node_set> seen;
    for(id_type i=0;i<labels.size();i++)
    {
        pair<unordered_map<id_type,node_set>::iterator,bool> ret = seen.insert(make_pair(labels[i],node_set()));
        ret.first->second.insert(i);
    }
    for(unordered_map<id_type,node_set>::iterator it = seen.begin(); it!= seen.end(); it++) communities.push_back(it->second);
}

void CDLib::convert_communities_to_labels(vector<node_set>& communities,vector<id_type>& labels)
{
    labels.clear();
    unordered_map<id_type,id_type> label_map;
    for(id_type i=0;i<communities.size();i++)
        for(node_set::const_iterator nit = communities[i].begin(); nit!= communities[i].end();nit++)
            label_map.insert(make_pair(*nit,i));
    labels.assign(label_map.size(),0);
    for(unordered_map<id_type,id_type>::iterator lmit = label_map.begin(); lmit != label_map.end();lmit++)
        labels[lmit->first] = lmit->second;
    
}

void CDLib::reindex_communities(const vector<id_type>& old_comms,vector<id_type>& new_comms){
    unordered_map<id_type,id_type> labelmap;
    new_comms.assign(old_comms.size(),0);
    id_type label_ctr = 0;
    for(id_type i=0;i<old_comms.size();i++){
        pair<unordered_map<id_type,id_type>::iterator,bool> ret = labelmap.insert(make_pair(old_comms[i],label_ctr));
        if(ret.second)label_ctr++;
        new_comms[i] = ret.first->second;
    }
}

bool CDLib::write_partition(const graph& g,const string& filepath,vector<node_set>& communities)
{
    ofstream ofs(filepath.c_str());
    if(!ofs.is_open()) return false;
    vector<id_type> labels(g.get_num_nodes(),0);
    for(id_type i=0;i<communities.size();i++)
        for(node_set::const_iterator nit = communities[i].begin(); nit!= communities[i].end();nit++)
            labels[*nit] = i;
    for(id_type i=0;i<labels.size();i++)
        ofs << g.get_node_label(i) << " " << labels[i] << endl;
    return true;
}

bool CDLib::write_partition_unlabelled(const graph& g,const string& filepath,vector<node_set>& communities)
{
    ofstream ofs(filepath.c_str());
    if(!ofs.is_open()) return false;
    vector<id_type> labels(g.get_num_nodes(),0);
    for(id_type i=0;i<communities.size();i++)
        for(node_set::const_iterator nit = communities[i].begin(); nit!= communities[i].end();nit++)
            labels[*nit] = i;
    for(id_type i=0;i<labels.size();i++)
        ofs << labels[i] << endl;
    return true;
}

bool CDLib::read_partition(const graph& g,const string& filepath,vector<node_set>& communities)
{
    vector<id_type> labels;
    labels.assign(g.get_num_nodes(),0);
    ifstream ifs(filepath);
    while(!ifs.eof())
    {
        string label; id_type comm_id;
        ifs >> label >> comm_id;
        if(label.size())labels[g.get_node_id(label)] = comm_id;
    }
    convert_labels_to_communities(labels,communities);
    return true;
}

id_type fix_communities(const graph& g1,const graph& g2,vector<node_set>& comms,vector<node_set>& comms_new,bool reverse)
{
    id_type common_elems=0;
    for(id_type i=0;i<comms.size();i++)
    {
        if(comms[i].size())
        {
            comms_new.push_back(node_set());
            for(node_set::iterator it=comms[i].begin();it != comms[i].end();it++)
            {
                id_type g1_node_id = *it;
                id_type g2_node_id = g2.get_node_id(g1.get_node_label(g1_node_id));
                if (g2_node_id < g2.get_num_nodes())
                {
                        comms_new[comms_new.size()-1].insert((reverse)?g2_node_id:g1_node_id);
                        common_elems++;
                }
            }
            if(comms_new[comms_new.size()-1].empty()) comms_new.pop_back();
        }
    }
    return common_elems;
}

void fix_communities_union(const graph& gc,node_set_string& unseen,unordered_map<string,id_type>& new_ids,vector<node_set>& comms,vector<node_set>& comms_new)
{
    for(id_type i=0;i<comms.size();i++)
    {
        if(comms[i].size())
        {
            comms_new.push_back(node_set());
            for(node_set::iterator it=comms[i].begin();it != comms[i].end();it++)
                comms_new[comms_new.size()-1].insert(new_ids[gc.get_node_label(*it)]);
        }
    }
    for(node_set_string::iterator it=unseen.begin();it != unseen.end();it++)
    {
        comms_new.push_back(node_set());
        comms_new[comms_new.size()-1].insert(new_ids[*it]);
    }
}



double CDLib::evolutionary_cluster_validation(const graph& g1, const graph& g2,vector<node_set>& comms1,vector<node_set>& comms2,double (*func)(id_type num_nodes, vector<node_set>& comms1, vector<node_set>& comms2))
{
    vector<node_set> comms_new1,comms_new2;
    id_type common1 = fix_communities(g1,g2,comms1,comms_new1,false);
    id_type common2 = fix_communities(g2,g1,comms2,comms_new2,true);
    if(common1 != common2) return 0;
    return func(common1,comms_new1,comms_new2);
}

double CDLib::evolutionary_cluster_validation_union(const graph& g1, const graph& g2,vector<node_set>& comms1,vector<node_set>& comms2,double (*func)(id_type num_nodes, vector<node_set>& comms1, vector<node_set>& comms2))
{
    vector<node_set> comms_new1,comms_new2;
    unordered_map<string,id_type> new_ids;
    node_set_string unseen1,unseen2;
    string g1_label,g2_label;
    for(id_type i=0; i<max(g1.get_num_nodes(),g2.get_num_nodes()) ;i++)
    {
        if(i<g1.get_num_nodes())
        {
            g1_label =  g1.get_node_label(i);
            if(g2.get_node_id(g1_label) >= g2.get_num_nodes()) unseen1.insert(g1_label);
            pair<unordered_map<string,id_type>::iterator,bool> ret = new_ids.insert(make_pair(g1_label,1));
            if(ret.second) ret.first->second = new_ids.size();
        }
        if(i<g2.get_num_nodes())
        {
            g2_label = g2.get_node_label(i);
            if(g1.get_node_id(g2_label) >= g1.get_num_nodes()) unseen2.insert(g2_label);
            pair<unordered_map<string,id_type>::iterator,bool> ret = new_ids.insert(make_pair(g2_label,1));
            if(ret.second) ret.first->second = new_ids.size();
        }
    }
    fix_communities_union(g2,unseen1,new_ids,comms1,comms_new1);
    fix_communities_union(g1,unseen2,new_ids,comms2,comms_new2);
    return func(new_ids.size(),comms_new1,comms_new2);
}

bool CDLib::is_member_of(id_type id, node_set& ns)
{
    return ns.find(id) != ns.end();
}

bool CDLib::in_same_comm(id_type i,id_type j,vector<node_set>& comms)
{
    for(id_type k=0;k<comms.size();k++)
        if(is_member_of(i,comms[k]) && is_member_of(j,comms[k])) 
            return true;
    return false;
}

id_type CDLib::in_comm(id_type i, vector<node_set>& comms)
{
    for(id_type k=0;k<comms.size();k++)
        if(is_member_of(i,comms[k]))
            return k;
    return comms.size();
}

double CDLib::degree_homogenity_test(graph& g, node_set & ns)
{
    graph sg(g.is_directed(),g.is_weighted());
    extract_subgraph(g,ns,sg);
    vector<id_type> degrees(sg.get_num_nodes(),0);
    for(id_type i=0;i<degrees.size();i++) degrees[i] = g.get_node_out_degree(i);
    statistics<id_type> degstats;
    compute_statistics(degrees,degstats);
    return static_cast<double>(degstats.variance)/static_cast<double>(degstats.mean_val);
}

double CDLib::kl_divergence(vector<double>& p,vector<double>& q)
{
    if(p.size() == q.size())
    {
        double kl_div = 0;
        for(id_type i=0;i<p.size();i++)
            if(p[i] && q[i]) 
                kl_div += p[i]*log(p[i]/q[i]);
        return kl_div;
    }
    return 1;
}

double CDLib::entropy_comparision_test(const graph& g,node_set& ns)
{
    graph subg(g.is_directed(),g.is_weighted());
    extract_subgraph(g,ns,subg);
    vector<double> qinit(subg.get_num_nodes(),1/static_cast<double>(subg.get_num_nodes())),outvec;
    id_type t = static_cast<id_type>(log(subg.get_num_nodes())/log(2));
    random_walk(subg,qinit,t,transform_func_row_stochastic,outvec);
    return kl_divergence(qinit,outvec);
}

void CDLib::get_community_graph(const graph&g, vector<node_set>& comms,graph& comm_graph)
{
    comm_graph.clear();
    vector<id_type> labels;
    convert_communities_to_labels(comms,labels);
    for(id_type i=0;i<comms.size();i++)comm_graph.add_node();   
    for(id_type i=0;i<comms.size();i++)
        for(node_set::iterator nit = comms[i].begin();nit!=comms[i].end();nit++)
            for(adjacent_edges_iterator aeit = g.out_edges_begin(*nit); aeit != g.out_edges_end(*nit); aeit++)
                    comm_graph.set_edge_weight(labels[*nit],labels[aeit->first],comm_graph.get_edge_weight(labels[*nit],labels[aeit->first]) + aeit->second);
}

void CDLib::compute_confusion_matrix_local(const graph& g,node_set& observed, node_set& truth, confusion_matrix_local& res)
{
    res.true_negatives = 0;
    res.true_positives = 0;
    res.false_negatives = 0;
    res.false_positives = 0;
    for(id_type i=0;i<g.get_num_nodes();i++)
    {
        if(is_member_of(i,observed) && is_member_of(i,truth)) res.true_positives++;
        else if (!is_member_of(i,observed) && !is_member_of(i,truth)) res.true_negatives++;
        else if (is_member_of(i,observed) && !is_member_of(i,truth)) res.false_positives++;
        else res.false_negatives++;
    }
    double truth_pos_prop = truth.size()/static_cast<double>(g.get_num_nodes());
    double truth_neg_prop = (g.get_num_nodes() - truth.size())/static_cast<double>(g.get_num_nodes());
    res.TPR = static_cast<double>(res.true_positives)/static_cast<double>(truth.size());
    res.TNR = static_cast<double>(res.true_negatives)/static_cast<double>(g.get_num_nodes() - truth.size());
    res.FPR = static_cast<double>(res.false_positives)/static_cast<double>(g.get_num_nodes() - truth.size());
    res.FNR = static_cast<double>(res.false_negatives)/static_cast<double>(truth.size());
    res.wFPR = truth_pos_prop*res.FPR + truth_neg_prop*res.FNR;
    res.wTPR = truth_pos_prop*res.TPR + truth_neg_prop*res.TNR;
    res.precision = static_cast<double>(res.true_positives)/static_cast<double>(observed.size());
    res.recall = res.TPR;
    res.f_score = (2*res.precision*res.recall)/(res.precision + res.recall);
}

void CDLib::compute_confusion_matrix_global(const graph&g, vector<node_set>& observed,vector<node_set>& truth,vector< vector<id_type> >& cmat)
{
    cmat.assign(observed.size(),vector<id_type>(truth.size(),0));
    for(id_type i=0;i<observed.size();i++)
        for(node_set::iterator nit = observed[i].begin(); nit != observed[i].end();nit++)
            for(id_type j=0;j<truth.size();j++)
                if(is_member_of(*nit,truth[j]))
                    cmat[i][j]++;
 
}