/* 
 * File:   lplabel_propagation.cpp
 * Author: bharath
 * 
 * Created on 3 May, 2012, 9:31 AM
 */


#include "label_propagation.h"
using namespace CDLib;

id_type CDLib::label_propagation_run(const graph&g,vector<id_type>& communities,id_type max_iters,lp_raghavan_2007* algoman)
{
    id_type num_iters=0;
    for(num_iters=0;num_iters<max_iters;num_iters++)
        if(algoman->do_iteration(g,num_iters)) break;
    algoman->finalize(g,num_iters,communities);
    return num_iters;
}

bool lp_raghavan_2007::is_node_lplabel_fixed(id_type node_id) const
{
    if(node_id > ids.size()) return false;
    else return (nodes_with_fixed_lplabels.find(node_id) != nodes_with_fixed_lplabels.end());
}

bool lp_raghavan_2007::fix_node_lplabel(id_type node_id)
{
    if(node_id > ids.size()) return false;
    nodes_with_fixed_lplabels.insert(node_id);
    return true;
}

bool lp_raghavan_2007::free_node_lplabel(id_type node_id)
{
    if(node_id > ids.size()) return false;
    nodes_with_fixed_lplabels.erase(node_id);
    return true;
}

bool lp_raghavan_2007::set_node_lplabel(id_type node_id,id_type lplabel)
{
    if(node_id > lplabels.size() && lplabel < lplabels.size() && is_node_lplabel_fixed(node_id)) return false;
    lplabels[node_id] = lplabel;
    return true;
}

id_type lp_raghavan_2007::get_node_lplabel(id_type node_id) const
{
    if(node_id > lplabels.size()) return lplabels.size();
    return lplabels[node_id];
}

void lp_raghavan_2007::post_node_assign(const graph& g, id_type current,id_type new_label,id_type num_iters, max_lplabel_container& max_labels) {}

void lp_raghavan_2007::post_iteration(const graph& g, id_type num_iters) {}

lp_raghavan_2007::lp_raghavan_2007(const graph& g,bool synchronous_val)
{
    synchronous = synchronous_val;
    ids.assign(g.get_num_nodes(),0);
    positions.assign(g.get_num_nodes(),0);
    lplabels.assign(g.get_num_nodes(),0);
    if(synchronous)lpnextlabels.assign(g.get_num_nodes(),0);
    for(id_type i=0;i<g.get_num_nodes();i++)
    {
        ids[i] = i;
        positions[i] = i;
        lplabels[i] = i;
        if(synchronous)lpnextlabels[i] = i;
    }
}

void lp_raghavan_2007::finalize(const graph& g,id_type num_iters,vector<id_type>& communities)
{
    communities.assign(g.get_num_nodes(),0);
    copy(lplabels.begin(),lplabels.end(),communities.begin());
}


bool lf_comp(const pair<id_type,double>& lhs,const pair<id_type,double>& rhs) { return lhs.second<rhs.second; }

void lp_raghavan_2007::get_max_lplabels(const graph& g,id_type current_node,max_lplabel_container& max_labels)
{
    lplabel_fitness_container lplabel_fitness_values;
    for(adjacent_edges_iterator aeit = g.in_edges_begin(current_node); aeit != g.in_edges_end(current_node); aeit++)
    {           
        double fitness_value = get_node_score(g,aeit->first)*get_label_fitness_for_edge(g,current_node,aeit->first,aeit->second);
        pair<lplabel_fitness_container::iterator,bool> ret = lplabel_fitness_values.insert(make_pair(lplabels[aeit->first],fitness_value));
        if(!ret.second) ret.first->second+= fitness_value;
    }
    lplabel_fitness_container::iterator mit = max_element(lplabel_fitness_values.begin(),lplabel_fitness_values.end(),lf_comp);
    double max_val = mit->second;
    while(mit != lplabel_fitness_values.end() && mit->second == max_val)
    {
        max_labels.push_back(mit->first);
        lplabel_fitness_values.erase(mit->first);
        mit = max_element(lplabel_fitness_values.begin(),lplabel_fitness_values.end(),lf_comp);
    }
}

bool lp_raghavan_2007::check_lplabels(const graph& g)
{
    for(id_type i=0; i<g.get_num_nodes(); i++)
    {
        vector<id_type> max_labels;
        get_max_lplabels(g,i,max_labels);
        if(find(max_labels.begin(),max_labels.end(),lplabels[i]) == max_labels.end())
            return false;
    }
    return true;
}

id_type lp_raghavan_2007::new_label_break_ties_randomly(const graph& g,id_type current_node,max_lplabel_container& max_labels)
{
   get_max_lplabels(g,current_node,max_labels);
   if(max_labels.empty()) return lplabels[current_node];
   CDLib::RandomGenerator<id_type> rnd_gen(0,max_labels.size()-1);
   return max_labels[rnd_gen.next()];
}

bool lp_raghavan_2007::do_iteration(const graph& g,id_type num_iters)
{
    for(vector<id_type>::iterator it=ids.begin(); it != ids.end();it++)
    {
        reorder(g,num_iters);
        max_lplabel_container max_labels;
        id_type new_label = new_label_break_ties_randomly(g,*it,max_labels);
        post_node_assign(g,*it,new_label,num_iters,max_labels);
        if(synchronous) lpnextlabels[*it] = new_label;
        else lplabels[*it] = new_label;
    }
    if(synchronous)
        for(id_type i=0;i<g.get_num_nodes();i++)
            lplabels[i] = lpnextlabels[i];
    post_iteration(g,num_iters);
    return check_lplabels(g);
}

void lp_raghavan_2007::reorder(const graph& g,id_type num_iters)
{
    random_shuffle(ids.begin(),ids.end());
    for(vector<id_type>::iterator it=ids.begin();it!=ids.end();it++)
        positions[*it] = it-ids.begin();
}

void lp_offensive_lpa::post_node_assign(const graph& g, id_type current,id_type new_label,id_type num_iters, max_lplabel_container& max_labels)
{
    lp_dyn_hop_att::post_node_assign(g,current,new_label,num_iters,max_labels);
    if(lplabels[current] != new_label)
        probabilities[current] = compute_probability_for_id(g,current,new_label,max_labels);
    
}

double lp_offensive_lpa::compute_probability_for_id(const graph& g,id_type current,id_type new_label,max_lplabel_container& max_labels)
{
    double prod_num=0;
    for(adjacent_edges_iterator aeit = g.in_edges_begin(current); aeit !=g.in_edges_end(current);aeit++)
    {
        if(lplabels[aeit->first] == new_label)prod_num += probabilities[aeit->first];
    }
    return prod_num/g.get_node_in_weight(current);
}

double lp_offensive_lpa::get_edge_weight_function(const graph&g,id_type from_id,id_type to_id,double edge_weight) { return edge_weight;}
double lp_offensive_lpa::get_label_fitness_for_edge(const graph& g, id_type from_id, id_type to_id,double edge_weight)
{
    return get_node_score(g,to_id)*(1-probabilities[to_id])*get_edge_weight_function(g,from_id,to_id,edge_weight);
}

double lp_offensive_lpa::get_node_probability(id_type id) const
{
    if(id > probabilities.size()) return 0.0;
    else return probabilities[id];
}

lp_offensive_lpa::lp_offensive_lpa(const graph& g,bool synchronous_val,double hop_att_max_val) : lp_dyn_hop_att(g,synchronous_val,hop_att_max_val)
{
    probabilities.assign(g.get_num_nodes(),(double)1/(double)g.get_num_nodes());
}

double lp_defensive_lpa::get_label_fitness_for_edge(const graph& g, id_type from_id, id_type to_id,double edge_weight)
{
    return get_node_score(g,to_id)*probabilities[to_id]*get_edge_weight_function(g,from_id,to_id,edge_weight);
}

double lp_defensive_lpa::compute_probability_for_id(const graph& g,id_type current,id_type new_label,max_lplabel_container& max_labels)
{
    double prob_denom=0,prod_num=0;
    for(adjacent_edges_iterator aeit = g.in_edges_begin(current); aeit !=g.in_edges_end(current);aeit++)
    {
        if(lplabels[aeit->first] == new_label)
        {
            prod_num += probabilities[aeit->first];
            prob_denom += aeit->second;
        }
    }
    return prod_num/prob_denom;
}

lp_defensive_lpa::lp_defensive_lpa(const graph& g,bool synchronous_val,double hop_att_max_val) : lp_offensive_lpa(g,synchronous_val,hop_att_max_val) {}

id_type CDLib::label_propagation_raghavan_2007(const graph& g,vector<id_type>& labels,id_type max_iters,bool synchronous)
{
    lp_raghavan_2007 algoman(g,synchronous);
    return label_propagation_run(g,labels,max_iters,&algoman);
}

id_type CDLib::label_propagation_leung_2009(const graph& g,vector<id_type>& labels,id_type max_iters,bool synchronous,double hop_att)
{
    lp_leung_2009 algoman(g,synchronous,hop_att);
    return label_propagation_run(g,labels,max_iters,&algoman);
}


id_type CDLib::label_propagation_dyn_hop_2010(const graph& g, vector<id_type>& labels,id_type max_iters,bool synchronous,double hop_att_max)
{
    lp_dyn_hop_att algoman(g,synchronous,hop_att_max);
    return label_propagation_run(g,labels,max_iters,&algoman);
}

id_type CDLib::label_propagation_olpa_2010(const graph& g, vector<id_type>&  labels,id_type max_iters,bool synchronous,double hop_att_max)
{
    lp_offensive_lpa algoman(g,synchronous,hop_att_max);
    return label_propagation_run(g,labels,max_iters,&algoman);
}

id_type CDLib::label_propagation_dlpa_2010(const graph& g, vector<id_type>&  labels,id_type max_iters,bool synchronous,double hop_att_max)
{
    lp_defensive_lpa algoman(g,synchronous,hop_att_max);
    return label_propagation_run(g,labels,max_iters,&algoman);
}

id_type CDLib::label_propagation_track_changes_2012(const graph& g, vector<id_type>&  labels,id_type max_iters,bool synchronous)
{
    lp_track_changes algoman(g,synchronous);
    return label_propagation_run(g,labels,max_iters,&algoman);
}