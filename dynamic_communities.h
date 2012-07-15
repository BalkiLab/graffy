/* 
 * File:   dynamic_communities.h
 * Author: bharath
 *
 * Created on 10 May, 2012, 11:24 AM
 */

#ifndef DYNAMIC_COMMUNITIES_H
#define	DYNAMIC_COMMUNITIES_H
#include "typedefs.h"
#include "profiler.h"
#include "graph.h"
#include "graphio.h"
#include "label_propagation.h"
namespace CDLib
{
    class evol_label_prop_new : public lp_track_changes
    {
    protected:
        double alpha;
        id_type num_iters_last;
        const vector<id_type>::iterator last_num_changes_it;
        virtual double get_label_fitness_for_edge(const graph& g, id_type from_id, id_type to_id,double edge_weight)
        {
            double history_fraction = 1.0-((*(last_num_changes_it + to_id))/num_iters_last);
            double snapshot_fraction = 1.0 - ((double)num_lplabel_changes[to_id]/(double)(num_iters_current+1));
            return (alpha*history_fraction + (1-alpha)*snapshot_fraction)*edge_weight;
        }
    public:
        evol_label_prop_new(const graph& g,bool synchronous_val,double alpha,id_type num_iters_last_val, vector<id_type>::iterator last_num_changes_begin) : lp_track_changes(g,synchronous_val),num_iters_last(num_iters_last_val),last_num_changes_it(last_num_changes_begin) {}
    };
    
    class evolutionary_label_propagation : public lp_dyn_hop_att
    {
    protected:
        bool activate_history_cost;
        bool copy_labels;
        double alpha;
        id_type num_iters_last;
        id_type num_iters_current;
        vector<id_type> last_num_lplabel_changes;
        evolutionary_label_propagation();
        virtual void post_iteration(const graph& g, id_type num_iters)
        {
            lp_dyn_hop_att::post_iteration(g,num_iters);
            num_iters_current = num_iters+1;
        }
        virtual double get_label_fitness_for_edge(const graph& g, id_type from_id, id_type to_id,double edge_weight)
        {
            if(activate_history_cost)
            {
                double history_fraction = 1.0-(last_num_lplabel_changes[to_id]/num_iters_last);
                double snapshot_fraction = 1.0 - ((double)num_lplabel_changes[to_id]/(double)(num_iters_current+1));
                return (alpha*history_fraction + (1-alpha)*snapshot_fraction)*edge_weight;
            }
            else return edge_weight;
        }
    public:
        evolutionary_label_propagation(const graph& g,bool synchronous_val,bool activate_history_cost_val,bool copy_labels_vals,double alpha_val,id_type num_iters_last_val,vector<id_type>& last_num_lplabel_changes_vec,vector<id_type>& last_labels_vec) 
            : lp_dyn_hop_att(g,synchronous_val,0),
            activate_history_cost(activate_history_cost_val),
            copy_labels(copy_labels_vals),
            alpha(alpha_val),
            num_iters_last(num_iters_last_val),
            num_iters_current(0),
            last_num_lplabel_changes(last_num_lplabel_changes_vec)
            {
                if(copy_labels)copy(last_labels_vec.begin(),last_labels_vec.end(),lplabels.begin());
            }
    };
    
    struct dynamic_lp_input
    {
        bool activate;
        bool copy_labels;
        double alpha;
        dynamic_lp_input(bool activate_val,bool copy_labels_val,double alpha_val) : activate(activate_val),copy_labels(copy_labels_val), alpha(alpha_val) {}
    };
    
    struct dynamic_lp_output
    {
       vector<graph> graphs;
       vector< vector<id_type> > lp_labels;
       vector<id_type> lp_iters;
       vector<double> lp_times;
       vector<double> book_times;
       vector< vector<id_type> > lp_num_label_changes;
    };
    
    double evolutionary_label_propagation_edgelists(const string& snapshot_filepath,bool directed, bool weighted,dynamic_lp_input& input, dynamic_lp_output& output);
    void new_evolutionary_label_propagation_edgelists(const string& snapshot_filepath,bool directed, bool weighted,double alpha,bool basic, dynamic_lp_output& output);
};

#endif	/* DYNAMIC_COMMUNITIES_H */

