/* 
 * File:   label_propagation.h
 * Author: bharath
 *
 * Created on 3 May, 2012, 9:31 AM
 */

#ifndef LABEL_PROPAGATION_H
#define	LABEL_PROPAGATION_H

#include "graph.h"
#include "community_tools.h"
using namespace std;

namespace CDLib
{
    typedef vector<id_type> max_lplabel_container;
    typedef unordered_map<id_type,double> lplabel_fitness_container;
    class lp_raghavan_2007
    {
    protected:
        bool synchronous;
        vector<id_type> positions;
        vector<id_type> ids;
        vector<id_type> lplabels;
        vector<id_type> lpnextlabels;
        unordered_set<id_type> nodes_with_fixed_lplabels;
        virtual double get_node_score(const graph& g, id_type node_id) { return 1.0; }
        virtual double get_label_fitness_for_edge(const graph& g, id_type from_id, id_type to_id,double edge_weight) { return edge_weight; }
        virtual void post_node_assign(const graph& g, id_type current,id_type new_label,id_type num_iters, max_lplabel_container& max_labels);
        virtual void post_iteration(const graph& g, id_type num_iters);
        virtual id_type new_label_break_ties_randomly(const graph& g,id_type current_node,max_lplabel_container& max_labels);
        virtual void get_max_lplabels(const graph& g,id_type current_node,max_lplabel_container& max_labels);
        virtual bool check_lplabels(const graph& g);
        virtual void reorder(const graph& g,id_type num_iters);
        lp_raghavan_2007();
    public:
        bool is_node_lplabel_fixed(id_type node_id) const;
        bool fix_node_lplabel(id_type node_id);
        bool free_node_lplabel(id_type node_id);
        bool set_node_lplabel(id_type node_id,id_type lplabel);
        id_type get_node_lplabel(id_type node_id) const;
        virtual bool do_iteration(const graph& g,id_type num_iters);
        virtual void finalize(const graph& g,id_type num_iters,vector<id_type>& communities);
        lp_raghavan_2007(const graph& g,bool synchronous_val);
    };
    
    id_type label_propagation_run(const graph&g,vector<id_type>& communities,id_type max_iters,lp_raghavan_2007* algoman);
    
    class lp_track_changes : public lp_raghavan_2007
    {
    protected:
        vector<id_type> num_lplabel_changes;
        id_type num_iters_current;
        virtual void post_iteration(const graph& g, id_type num_iters)
        {
            num_iters_current = num_iters+1;
        }
        virtual void post_node_assign(const graph& g, id_type current,id_type new_label,id_type num_iters, max_lplabel_container& max_labels)
        {
            if(lplabels[current] != new_label)
                num_lplabel_changes[current]++;
        }
        virtual double get_label_fitness_for_edge(const graph& g, id_type from_id, id_type to_id,double edge_weight)
        {
            double snapshot_fraction = 1.0 - ((double)num_lplabel_changes[to_id]/(double)(num_iters_current+1));
            return snapshot_fraction*edge_weight;
        }
    public:
        vector<id_type>::const_iterator label_changes_begin() const { return num_lplabel_changes.begin();}
        vector<id_type>::const_iterator label_changes_end() const { return num_lplabel_changes.end();}
        id_type get_node_num_lplabel_changes(id_type node_id) const
        {
            if(node_id > num_lplabel_changes.size()) return 0;
                    return num_lplabel_changes[node_id];
        }
        lp_track_changes(const graph& g,bool synchronous_val): lp_raghavan_2007(g,synchronous_val), num_lplabel_changes(g.get_num_nodes(),0), num_iters_current(0) {}
    };
    
    
    class lp_leung_2009 : public lp_raghavan_2007
    {
    protected:
        double hop_att;
        vector<double> scores;
        virtual double get_node_score(const graph& g, id_type node_id) {  return (scores[node_id] - hop_att); }
        virtual void post_node_assign(const graph& g, id_type current,id_type new_label,id_type num_iters, max_lplabel_container& max_labels)
        {
            lp_raghavan_2007::post_node_assign(g,current,new_label,num_iters,max_labels);
            if(lplabels[current] != new_label)
                scores[current] = compute_score_for_id(g,current,new_label,max_labels);
        }
        virtual double compute_score_for_id(const graph& g,id_type current_node,id_type label,max_lplabel_container& max_lplabels)
        {
            double ret_val = -numeric_limits<double>::infinity();
            for(adjacent_edges_iterator aeit = g.in_edges_begin(current_node); aeit != g.in_edges_end(current_node); aeit++)
                if(lplabels[aeit->first] == label && ret_val < scores[aeit->first])
                    ret_val = scores[aeit->first];
            return ret_val;
        }
    public:
        double get_node_score_value(id_type node_id) const
        {
            if(node_id > scores.size()) return 0.0;
                    return scores[node_id];
        }
        bool set_node_score_value(id_type node_id,double score)
        {
            if (node_id > scores.size()) return false;
            scores[node_id] = score;
            return true;
        }
        lp_leung_2009(const graph& g,bool synchronous_val,double hop_att_val) :
                lp_raghavan_2007(g,synchronous_val),
                hop_att(hop_att_val),
                scores(vector<double>(g.get_num_nodes(),0)) {}

    };
   

    class lp_dyn_hop_att : public lp_track_changes
    {
    protected:
        double hop_att_max;
        double current_hop_att;
        vector<double> distances;
        virtual void post_iteration(const graph& g, id_type num_iters)
        {
            lp_track_changes::post_iteration(g,num_iters);
            double prop_change = 0;
            for(id_type i=0;i<num_lplabel_changes.size();i++)
                prop_change+=num_lplabel_changes[i];
            prop_change /= g.get_num_nodes();
            current_hop_att = (prop_change > hop_att_max) ? 0 : current_hop_att;
        }
        virtual void post_node_assign(const graph& g, id_type current,id_type new_label,id_type num_iters, max_lplabel_container& max_labels)
        {
            lp_track_changes::post_node_assign(g,current,new_label,num_iters,max_labels);
            if(lplabels[current] != new_label)
            {
                distances[current] = compute_distance_for_id(g,current,new_label,max_labels);
            }
        }
        virtual double compute_distance_for_id(const graph& g,id_type current_node,id_type new_label,max_lplabel_container& max_lplabels)
        {
            double ret_val = numeric_limits<double>::infinity();
            for(adjacent_edges_iterator aeit = g.in_edges_begin(current_node); aeit != g.in_edges_end(current_node); aeit++)
                if(lplabels[aeit->first] == new_label && ret_val > distances[aeit->first])
                    ret_val = distances[aeit->first];
            return ret_val + 1;
        }
        virtual double get_node_score(const graph& g, id_type node_id) { return 1 - (current_hop_att * distances[node_id]); }
        lp_dyn_hop_att();
    public:
        double get_node_distance(id_type node_id) const
        {
            if(node_id > distances.size()) return 0.0;
            return distances[node_id];
        }
        bool set_node_distance(id_type node_id,double distance)
        {
            if(node_id > distances.size()) return false;
            distances[node_id] = distance;
            return true;
        }
        lp_dyn_hop_att(const graph& g,bool synchronous_val,double hop_att_max_val) : 
                        lp_track_changes(g,synchronous), 
                        hop_att_max(hop_att_max_val), 
                        current_hop_att(0),
                        distances(vector<double>(g.get_num_nodes(),0)) {}
    };
    
    class lp_offensive_lpa : public lp_dyn_hop_att
    {
    protected:
        lp_offensive_lpa();
        vector<double> probabilities;
        virtual void post_node_assign(const graph& g, id_type current,id_type new_label,id_type num_iters, max_lplabel_container& max_labels);
        virtual double compute_probability_for_id(const graph& g,id_type current,id_type new_label,max_lplabel_container& max_labels);
        virtual double get_edge_weight_function(const graph&g,id_type from_id,id_type to_id,double edge_weight);
        virtual double get_label_fitness_for_edge(const graph& g, id_type from_id, id_type to_id,double edge_weight);
        
    public:
        double get_node_probability(id_type id) const;
        bool set_node_probability(id_type id,double probability);
        lp_offensive_lpa(const graph& g,bool synchronous_val,double hop_att_max_val);
    };
    
    class lp_defensive_lpa : public lp_offensive_lpa
    {
    protected:
        lp_defensive_lpa();
        virtual double compute_probability_for_id(const graph& g,id_type current,id_type new_label,max_lplabel_container& max_labels);
        virtual double get_label_fitness_for_edge(const graph& g, id_type from_id, id_type to_id,double edge_weight);
    public:
        lp_defensive_lpa(const graph& g,bool synchronous_val,double hop_att_max_val);
    };
    
    template<class node_preference_function,class edge_weight_function>
    class extendable_lp_raghavan_2007 : public lp_raghavan_2007
    {
    protected:
        node_preference_function node_pref_func;
        edge_weight_function edge_weight_func;
        virtual double get_label_fitness_for_edge(const graph& g, id_type from_id, id_type to_id,double edge_weight)
        {
            return get_node_score(g,from_id)*node_pref_func(g,to_id)*edge_weight_func(g,from_id,to_id,edge_weight);
        }
        extendable_lp_raghavan_2007();
    public:
        extendable_lp_raghavan_2007(const graph& g, bool synchronous_val,node_preference_function& node_pref_f, edge_weight_function& edge_wt_f) : lp_raghavan_2007(g,synchronous_val),node_pref_func(node_pref_f),edge_weight_func(edge_wt_f) {}
    };
    
    
    template<class node_preference_function,class edge_weight_function>
    class extendable_lp_leung_2009 : public lp_leung_2009
    {
    protected:
        node_preference_function node_pref_func;
        edge_weight_function edge_weight_func;
        virtual double get_label_fitness_for_edge(const graph& g, id_type from_id, id_type to_id,double edge_weight)
        {
            return get_node_score(g,from_id)*node_pref_func(g,to_id)*edge_weight_func(g,from_id,to_id,edge_weight);
        }
        extendable_lp_leung_2009();
    public:
        extendable_lp_leung_2009(const graph& g, bool synchronous_val,double hop_att_val,node_preference_function& node_pref_f, edge_weight_function& edge_wt_f) : lp_leung_2009(g,synchronous,hop_att_val),node_pref_func(node_pref_f),edge_weight_func(edge_wt_f) {}
    };
   
    
    template<class node_preference_function_object,class edge_weight_function_object>
    class extendable_lp_dyn_hop_att : public lp_dyn_hop_att
    {
    protected:
        node_preference_function_object node_pref_func;
        edge_weight_function_object edge_weight_func;
        virtual double get_label_fitness_for_edge(const graph& g, id_type from_id, id_type to_id,double edge_weight)
        {
            return get_node_score(g,from_id)*node_pref_func(g,to_id)*edge_weight_func(g,from_id,to_id,edge_weight);
        }
        extendable_lp_dyn_hop_att();
    public:
        extendable_lp_dyn_hop_att(const graph& g, bool synchronous_val,double hop_att_max_val,node_preference_function_object& node_pref_f, edge_weight_function_object& edge_wt_f) : lp_dyn_hop_att(g,synchronous,hop_att_max_val),node_pref_func(node_pref_f),edge_weight_func(edge_wt_f) {}
    };
    
    template<class edge_weight_function_object>
    class extendable_lp_offensive_lpa : public lp_offensive_lpa
    {
    protected:
        edge_weight_function_object edge_weight_func;
        virtual double get_edge_weight_function(const graph& g, id_type from_id, id_type to_id,double edge_weight)
        {
            return edge_weight_func(g,from_id,to_id,edge_weight);
        }
        extendable_lp_offensive_lpa();
    public:
        extendable_lp_offensive_lpa(const graph& g, bool synchronous_val,double hop_att_max_val,edge_weight_function_object& edge_wt_f) : lp_offensive_lpa(g,synchronous,hop_att_max_val),edge_weight_func(edge_wt_f) {}
    };
    
    template<class edge_weight_function_object>
    class extendable_lp_defensive_lpa : public lp_defensive_lpa
    {
    protected:
        edge_weight_function_object edge_weight_func;
        virtual double get_edge_weight_function(const graph& g, id_type from_id, id_type to_id,double edge_weight)
        {
            return edge_weight_func(g,from_id,to_id,edge_weight);
        }
        extendable_lp_defensive_lpa();
    public:
        extendable_lp_defensive_lpa(const graph& g, bool synchronous_val,double hop_att_max_val,edge_weight_function_object& edge_wt_f) : lp_defensive_lpa(g,synchronous,hop_att_max_val),edge_weight_func(edge_wt_f) {}
    };
    
    id_type label_propagation_raghavan_2007(const graph& g, vector<id_type>& labels,id_type max_iters,bool synchronous);
    id_type label_propagation_leung_2009(const graph& g, vector<id_type>&  labels,id_type max_iters,bool synchronous,double hop_att);
    id_type label_propagation_dyn_hop_2010(const graph& g, vector<id_type>& labels,id_type max_iters,bool synchronous,double hop_att_max);
    id_type label_propagation_olpa_2010(const graph& g, vector<id_type>&  labels,id_type max_iters,bool synchronous,double hop_att_max);
    id_type label_propagation_dlpa_2010(const graph& g, vector<id_type>&  labels,id_type max_iters,bool synchronous,double hop_att_max);
    id_type label_propagation_track_changes_2012(const graph& g, vector<id_type>&  labels,id_type max_iters,bool synchronous);
   
};

#endif	/* lplabel_PROPAGATION_H */

