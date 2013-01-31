/* 
 * File:   community.h
 * Author: bharath
 *
 * Created on 31 January, 2013, 12:24 AM
 */

#ifndef COMMUNITY_H
#define	COMMUNITY_H

#include "graph.h"
#include "community_tools.h"
#include "graphio.h"
#include "graph_operations.h"
#include "profiler.h"
#include "utility.h"
#include "centrality.h"
#include "paths_and_components.h"

namespace CDLib {

    bool local_community_clauset(const graph& g, id_type src, size_t k, node_set& output);
    bool local_community_clauset_modified(const graph& g, id_type src, size_t k, node_set& output);
    bool LWP_2006(const graph& g, id_type src, node_set& output);
    void VD_2011(const graph& g, id_type src, node_set& output);
    bool Bagrow_2007(const graph& g, id_type src, node_set& output);
    bool My_Algorithm(const graph& g, id_type src, node_set& output);
    bool CZR(const graph& g, id_type src, node_set& output);
    bool CZR_Beta(const graph& g, id_type src, node_set& output);

    typedef vector< vector<node_set> > dendrogram;
    void girvan_newman_2002(const graph& g, dendrogram& dendro);
    void radicchi_et_al_2004(const graph& g, dendrogram & dendro);

    typedef vector<id_type> max_lplabel_container;
    typedef unordered_map<id_type, double> lplabel_fitness_container;

    class lp_raghavan_2007 {
    protected:
        bool synchronous;
        vector<id_type> positions;
        vector<id_type> ids;
        vector<id_type> lplabels;
        vector<id_type> lpnextlabels;
        unordered_set<id_type> nodes_with_fixed_lplabels;

        virtual double get_node_score(const graph& g, id_type node_id) {
            return 1.0;
        }

        virtual double get_label_fitness_for_edge(const graph& g, id_type from_id, id_type to_id, double edge_weight) {
            return edge_weight;
        }
        virtual void post_node_assign(const graph& g, id_type current, id_type new_label, id_type num_iters, max_lplabel_container& max_labels);
        virtual void post_iteration(const graph& g, id_type num_iters);
        virtual id_type new_label_break_ties_randomly(const graph& g, id_type current_node, max_lplabel_container& max_labels);
        virtual void get_max_lplabels(const graph& g, id_type current_node, max_lplabel_container& max_labels);
        virtual bool check_lplabels(const graph& g);
        virtual void reorder(const graph& g, id_type num_iters);
        lp_raghavan_2007();
    public:
        bool is_node_lplabel_fixed(id_type node_id) const;
        bool fix_node_lplabel(id_type node_id);
        bool free_node_lplabel(id_type node_id);
        bool set_node_lplabel(id_type node_id, id_type lplabel);
        id_type get_node_lplabel(id_type node_id) const;
        virtual bool do_iteration(const graph& g, id_type num_iters);
        virtual void finalize(const graph& g, id_type num_iters, vector<id_type>& communities);
        lp_raghavan_2007(const graph& g, bool synchronous_val);
    };

    id_type label_propagation_run(const graph&g, vector<id_type>& communities, id_type max_iters, lp_raghavan_2007* algoman);

    class lp_track_changes : public lp_raghavan_2007 {
    protected:
        vector<id_type> num_lplabel_changes;
        id_type num_iters_current;

        virtual void post_iteration(const graph& g, id_type num_iters) {
            num_iters_current = num_iters + 1;
        }

        virtual void post_node_assign(const graph& g, id_type current, id_type new_label, id_type num_iters, max_lplabel_container& max_labels) {
            if (lplabels[current] != new_label)
                num_lplabel_changes[current]++;
        }

        virtual double get_label_fitness_for_edge(const graph& g, id_type from_id, id_type to_id, double edge_weight) {
            double snapshot_fraction = 1.0 - ((double) num_lplabel_changes[to_id] / (double) (num_iters_current + 1));
            return snapshot_fraction*edge_weight;
        }
    public:

        vector<id_type>::const_iterator label_changes_begin() const {
            return num_lplabel_changes.begin();
        }

        vector<id_type>::const_iterator label_changes_end() const {
            return num_lplabel_changes.end();
        }

        id_type get_node_num_lplabel_changes(id_type node_id) const {
            if (node_id > num_lplabel_changes.size()) return 0;
            return num_lplabel_changes[node_id];
        }

        lp_track_changes(const graph& g, bool synchronous_val) : lp_raghavan_2007(g, synchronous_val), num_lplabel_changes(g.get_num_nodes(), 0), num_iters_current(0) {
        }
    };

    class lp_leung_2009 : public lp_raghavan_2007 {
    protected:
        double hop_att;
        vector<double> scores;

        virtual double get_node_score(const graph& g, id_type node_id) {
            return (scores[node_id] - hop_att);
        }

        virtual void post_node_assign(const graph& g, id_type current, id_type new_label, id_type num_iters, max_lplabel_container& max_labels) {
            lp_raghavan_2007::post_node_assign(g, current, new_label, num_iters, max_labels);
            if (lplabels[current] != new_label)
                scores[current] = compute_score_for_id(g, current, new_label, max_labels);
        }

        virtual double compute_score_for_id(const graph& g, id_type current_node, id_type label, max_lplabel_container& max_lplabels) {
            double ret_val = -numeric_limits<double>::infinity();
            for (adjacent_edges_iterator aeit = g.in_edges_begin(current_node); aeit != g.in_edges_end(current_node); aeit++)
                if (lplabels[aeit->first] == label && ret_val < scores[aeit->first])
                    ret_val = scores[aeit->first];
            return ret_val;
        }
    public:

        double get_node_score_value(id_type node_id) const {
            if (node_id > scores.size()) return 0.0;
            return scores[node_id];
        }

        bool set_node_score_value(id_type node_id, double score) {
            if (node_id > scores.size()) return false;
            scores[node_id] = score;
            return true;
        }

        lp_leung_2009(const graph& g, bool synchronous_val, double hop_att_val) :
        lp_raghavan_2007(g, synchronous_val),
        hop_att(hop_att_val),
        scores(vector<double>(g.get_num_nodes(), 0)) {
        }

    };

    class lp_dyn_hop_att : public lp_track_changes {
    protected:
        double hop_att_max;
        double current_hop_att;
        vector<double> distances;

        virtual void post_iteration(const graph& g, id_type num_iters) {
            lp_track_changes::post_iteration(g, num_iters);
            double prop_change = 0;
            for (id_type i = 0; i < num_lplabel_changes.size(); i++)
                prop_change += num_lplabel_changes[i];
            prop_change /= g.get_num_nodes();
            current_hop_att = (prop_change > hop_att_max) ? 0 : current_hop_att;
        }

        virtual void post_node_assign(const graph& g, id_type current, id_type new_label, id_type num_iters, max_lplabel_container& max_labels) {
            lp_track_changes::post_node_assign(g, current, new_label, num_iters, max_labels);
            if (lplabels[current] != new_label) {
                distances[current] = compute_distance_for_id(g, current, new_label, max_labels);
            }
        }

        virtual double compute_distance_for_id(const graph& g, id_type current_node, id_type new_label, max_lplabel_container& max_lplabels) {
            double ret_val = numeric_limits<double>::infinity();
            for (adjacent_edges_iterator aeit = g.in_edges_begin(current_node); aeit != g.in_edges_end(current_node); aeit++)
                if (lplabels[aeit->first] == new_label && ret_val > distances[aeit->first])
                    ret_val = distances[aeit->first];
            return ret_val + 1;
        }

        virtual double get_node_score(const graph& g, id_type node_id) {
            return 1 - (current_hop_att * distances[node_id]);
        }
        lp_dyn_hop_att();
    public:

        double get_node_distance(id_type node_id) const {
            if (node_id > distances.size()) return 0.0;
            return distances[node_id];
        }

        bool set_node_distance(id_type node_id, double distance) {
            if (node_id > distances.size()) return false;
            distances[node_id] = distance;
            return true;
        }

        lp_dyn_hop_att(const graph& g, bool synchronous_val, double hop_att_max_val) :
        lp_track_changes(g, synchronous),
        hop_att_max(hop_att_max_val),
        current_hop_att(0),
        distances(vector<double>(g.get_num_nodes(), 0)) {
        }
    };

    class lp_offensive_lpa : public lp_dyn_hop_att {
    protected:
        lp_offensive_lpa();
        vector<double> probabilities;
        virtual void post_node_assign(const graph& g, id_type current, id_type new_label, id_type num_iters, max_lplabel_container& max_labels);
        virtual double compute_probability_for_id(const graph& g, id_type current, id_type new_label, max_lplabel_container& max_labels);
        virtual double get_edge_weight_function(const graph&g, id_type from_id, id_type to_id, double edge_weight);
        virtual double get_label_fitness_for_edge(const graph& g, id_type from_id, id_type to_id, double edge_weight);

    public:
        double get_node_probability(id_type id) const;
        bool set_node_probability(id_type id, double probability);
        lp_offensive_lpa(const graph& g, bool synchronous_val, double hop_att_max_val);
    };

    class lp_defensive_lpa : public lp_offensive_lpa {
    protected:
        lp_defensive_lpa();
        virtual double compute_probability_for_id(const graph& g, id_type current, id_type new_label, max_lplabel_container& max_labels);
        virtual double get_label_fitness_for_edge(const graph& g, id_type from_id, id_type to_id, double edge_weight);
    public:
        lp_defensive_lpa(const graph& g, bool synchronous_val, double hop_att_max_val);
    };

    template<class node_preference_function, class edge_weight_function>
    class extendable_lp_raghavan_2007 : public lp_raghavan_2007 {
    protected:
        node_preference_function node_pref_func;
        edge_weight_function edge_weight_func;

        virtual double get_label_fitness_for_edge(const graph& g, id_type from_id, id_type to_id, double edge_weight) {
            return get_node_score(g, from_id) * node_pref_func(g, to_id) * edge_weight_func(g, from_id, to_id, edge_weight);
        }
        extendable_lp_raghavan_2007();
    public:

        extendable_lp_raghavan_2007(const graph& g, bool synchronous_val, node_preference_function& node_pref_f, edge_weight_function& edge_wt_f) : lp_raghavan_2007(g, synchronous_val), node_pref_func(node_pref_f), edge_weight_func(edge_wt_f) {
        }
    };

    template<class node_preference_function, class edge_weight_function>
    class extendable_lp_leung_2009 : public lp_leung_2009 {
    protected:
        node_preference_function node_pref_func;
        edge_weight_function edge_weight_func;

        virtual double get_label_fitness_for_edge(const graph& g, id_type from_id, id_type to_id, double edge_weight) {
            return get_node_score(g, from_id) * node_pref_func(g, to_id) * edge_weight_func(g, from_id, to_id, edge_weight);
        }
        extendable_lp_leung_2009();
    public:

        extendable_lp_leung_2009(const graph& g, bool synchronous_val, double hop_att_val, node_preference_function& node_pref_f, edge_weight_function& edge_wt_f) : lp_leung_2009(g, synchronous, hop_att_val), node_pref_func(node_pref_f), edge_weight_func(edge_wt_f) {
        }
    };

    template<class node_preference_function_object, class edge_weight_function_object>
    class extendable_lp_dyn_hop_att : public lp_dyn_hop_att {
    protected:
        node_preference_function_object node_pref_func;
        edge_weight_function_object edge_weight_func;

        virtual double get_label_fitness_for_edge(const graph& g, id_type from_id, id_type to_id, double edge_weight) {
            return get_node_score(g, from_id) * node_pref_func(g, to_id) * edge_weight_func(g, from_id, to_id, edge_weight);
        }
        extendable_lp_dyn_hop_att();
    public:

        extendable_lp_dyn_hop_att(const graph& g, bool synchronous_val, double hop_att_max_val, node_preference_function_object& node_pref_f, edge_weight_function_object& edge_wt_f) : lp_dyn_hop_att(g, synchronous, hop_att_max_val), node_pref_func(node_pref_f), edge_weight_func(edge_wt_f) {
        }
    };

    template<class edge_weight_function_object>
    class extendable_lp_offensive_lpa : public lp_offensive_lpa {
    protected:
        edge_weight_function_object edge_weight_func;

        virtual double get_edge_weight_function(const graph& g, id_type from_id, id_type to_id, double edge_weight) {
            return edge_weight_func(g, from_id, to_id, edge_weight);
        }
        extendable_lp_offensive_lpa();
    public:

        extendable_lp_offensive_lpa(const graph& g, bool synchronous_val, double hop_att_max_val, edge_weight_function_object& edge_wt_f) : lp_offensive_lpa(g, synchronous, hop_att_max_val), edge_weight_func(edge_wt_f) {
        }
    };

    template<class edge_weight_function_object>
    class extendable_lp_defensive_lpa : public lp_defensive_lpa {
    protected:
        edge_weight_function_object edge_weight_func;

        virtual double get_edge_weight_function(const graph& g, id_type from_id, id_type to_id, double edge_weight) {
            return edge_weight_func(g, from_id, to_id, edge_weight);
        }
        extendable_lp_defensive_lpa();
    public:

        extendable_lp_defensive_lpa(const graph& g, bool synchronous_val, double hop_att_max_val, edge_weight_function_object& edge_wt_f) : lp_defensive_lpa(g, synchronous, hop_att_max_val), edge_weight_func(edge_wt_f) {
        }
    };

    id_type label_propagation_raghavan_2007(const graph& g, vector<id_type>& labels, id_type max_iters, bool synchronous);
    id_type label_propagation_leung_2009(const graph& g, vector<id_type>& labels, id_type max_iters, bool synchronous, double hop_att);
    id_type label_propagation_dyn_hop_2010(const graph& g, vector<id_type>& labels, id_type max_iters, bool synchronous, double hop_att_max);
    id_type label_propagation_olpa_2010(const graph& g, vector<id_type>& labels, id_type max_iters, bool synchronous, double hop_att_max);
    id_type label_propagation_dlpa_2010(const graph& g, vector<id_type>& labels, id_type max_iters, bool synchronous, double hop_att_max);
    id_type label_propagation_track_changes_2012(const graph& g, vector<id_type>& labels, id_type max_iters, bool synchronous);

    class evol_label_prop_new : public lp_track_changes {
    protected:
        double alpha;
        id_type num_iters_last;
        const vector<id_type>::iterator last_num_changes_it;

        virtual double get_label_fitness_for_edge(const graph& g, id_type from_id, id_type to_id, double edge_weight) {
            double history_fraction = 1.0 - ((*(last_num_changes_it + to_id)) / num_iters_last);
            double snapshot_fraction = 1.0 - ((double) num_lplabel_changes[to_id] / (double) (num_iters_current + 1));
            return (alpha * history_fraction + (1 - alpha) * snapshot_fraction)*edge_weight;
        }
    public:

        evol_label_prop_new(const graph& g, bool synchronous_val, double alpha, id_type num_iters_last_val, vector<id_type>::iterator last_num_changes_begin) : lp_track_changes(g, synchronous_val), num_iters_last(num_iters_last_val), last_num_changes_it(last_num_changes_begin) {
        }
    };

    class evolutionary_label_propagation : public lp_dyn_hop_att {
    protected:
        bool activate_history_cost;
        bool copy_labels;
        double alpha;
        id_type num_iters_last;
        id_type num_iters_current;
        vector<id_type> last_num_lplabel_changes;
        evolutionary_label_propagation();

        virtual void post_iteration(const graph& g, id_type num_iters) {
            lp_dyn_hop_att::post_iteration(g, num_iters);
            num_iters_current = num_iters + 1;
        }

        virtual double get_label_fitness_for_edge(const graph& g, id_type from_id, id_type to_id, double edge_weight) {
            if (activate_history_cost) {
                double history_fraction = 1.0 - (last_num_lplabel_changes[to_id] / num_iters_last);
                double snapshot_fraction = 1.0 - ((double) num_lplabel_changes[to_id] / (double) (num_iters_current + 1));
                return (alpha * history_fraction + (1 - alpha) * snapshot_fraction)*edge_weight;
            } else return edge_weight;
        }
    public:

        evolutionary_label_propagation(const graph& g, bool synchronous_val, bool activate_history_cost_val, bool copy_labels_vals, double alpha_val, id_type num_iters_last_val, vector<id_type>& last_num_lplabel_changes_vec, vector<id_type>& last_labels_vec)
        : lp_dyn_hop_att(g, synchronous_val, 0),
        activate_history_cost(activate_history_cost_val),
        copy_labels(copy_labels_vals),
        alpha(alpha_val),
        num_iters_last(num_iters_last_val),
        num_iters_current(0),
        last_num_lplabel_changes(last_num_lplabel_changes_vec) {
            if (copy_labels)copy(last_labels_vec.begin(), last_labels_vec.end(), lplabels.begin());
        }
    };

    struct dynamic_lp_input {
        bool activate;
        bool copy_labels;
        double alpha;

        dynamic_lp_input(bool activate_val, bool copy_labels_val, double alpha_val) : activate(activate_val), copy_labels(copy_labels_val), alpha(alpha_val) {
        }
    };

    struct dynamic_lp_output {
        vector<graph> graphs;
        vector< vector<id_type> > lp_labels;
        vector<id_type> lp_iters;
        vector<double> lp_times;
        vector<double> book_times;
        vector< vector<id_type> > lp_num_label_changes;
    };

    double evolutionary_label_propagation_edgelists(const string& snapshot_filepath, bool directed, bool weighted, dynamic_lp_input& input, dynamic_lp_output& output);
    void new_evolutionary_label_propagation_edgelists(const string& snapshot_filepath, bool directed, bool weighted, double alpha, bool basic, dynamic_lp_output& output);

    class bgll_objective {
    public:
        virtual void init(const graph& orig_graph, const vector<vector<id_type> >& hiercomms, const graph& curr_graph, const vector<id_type>& curr_comms) = 0;
        virtual double objval(const graph& curr_graph, const vector<id_type>& curr_comms) const = 0;
        virtual double compute_gain(const graph& curr_graph, const vector<id_type>& labels, id_type vertex, id_type dst_comm) const = 0;
        virtual void detach_node(const graph& curr_graph, const vector<id_type>& labels, id_type vertex) = 0;
        virtual void attach_node(const graph& curr_graph, const vector<id_type>& labels, id_type vertex, id_type dst_comm) = 0;
    };

    class bgll_objective_w_internal : public bgll_objective {
    protected:
        unordered_map<id_type, double> internal_weight;
        unordered_map<id_type, double> node_link_weights;
        double resolution_param;
        double self_weight;
        double weight_inside;

        inline void init_internal_params() {
            internal_weight.clear();
            node_link_weights.clear();
            self_weight = 0;
            weight_inside = 0;
        }

        inline void fill_internal_edges_for_vertex(const graph& curr_graph, const vector<id_type>& curr_comms, id_type vertex) {
            for (adjacent_edges_iterator aeit = curr_graph.out_edges_begin(vertex); aeit != curr_graph.out_edges_end(vertex); aeit++)
                if (curr_comms[vertex] == curr_comms[aeit->first])
                    map_insert_and_increment<id_type,double>(internal_weight, curr_comms[vertex], aeit->second);
        }

        inline void process_vertex_links(const graph& curr_graph, const vector<id_type>& labels, id_type vertex) {
            node_link_weights.clear();
            self_weight = 0;
            weight_inside = 0;
            for (adjacent_edges_iterator aeit = curr_graph.out_edges_begin(vertex); aeit != curr_graph.out_edges_end(vertex); aeit++) {
                if (labels[aeit->first] != labels[vertex]) {
                    pair < unordered_map<id_type, double>::iterator, bool> ret = node_link_weights.insert(make_pair(labels[aeit->first], 0));
                    ret.first->second += aeit->second;
                } else if (vertex != aeit->first) {
                    weight_inside += aeit->second;
                } else self_weight = aeit->second;
            }
        }

        inline void update_internal_edges_on_vertex_detach(const graph& curr_graph, const vector<id_type>& labels, id_type vertex) {
            process_vertex_links(curr_graph, labels, vertex);
            unordered_map<id_type, double>::iterator it = internal_weight.find(labels[vertex]);
            if (it != internal_weight.end()) it->second -= (2 * weight_inside + self_weight);
        }

        inline void update_internal_edges_on_vertex_attach(id_type dst_comm) {
            unordered_map<id_type, double>::iterator it_in = internal_weight.find(dst_comm), it_nl = node_link_weights.find(dst_comm);
            double int_wt = self_weight;
            if (it_nl != node_link_weights.end()) int_wt += (2 * it_nl->second);
            if (it_in != internal_weight.end()) it_in->second += int_wt;
            else internal_weight.insert(make_pair(dst_comm, int_wt));
        }

        inline double objval_internal_term_running(const graph& curr_graph, id_type comm_id) const {
            double int_weight = 0;
            unordered_map<id_type, double>::const_iterator cit = internal_weight.find(comm_id);
            if (cit != internal_weight.end()) int_weight = cit->second;
            return ((resolution_param * int_weight) / (2 * curr_graph.get_total_weight()));
        }

        inline double gain_internal_term(id_type dst_comm) const {
            unordered_map<id_type, double>::const_iterator it = node_link_weights.find(dst_comm);
            if (it != node_link_weights.end()) return (resolution_param * ((2 * it->second) + self_weight));
            else return 0;
        }
    public:
        bgll_objective_w_internal() : resolution_param(1) {
        }

        bgll_objective_w_internal(double resolution_p) : resolution_param(resolution_p) {
        }
        virtual void init(const graph& orig_graph, const vector<vector<id_type> >& hiercomms, const graph& curr_graph, const vector<id_type>& curr_comms) = 0;
        virtual double objval(const graph& curr_graph, const vector<id_type>& curr_comms) const = 0;
        virtual double compute_gain(const graph& curr_graph, const vector<id_type>& labels, id_type vertex, id_type dst_comm) const = 0;
        virtual void detach_node(const graph& curr_graph, const vector<id_type>& labels, id_type vertex) = 0;
        virtual void attach_node(const graph& curr_graph, const vector<id_type>& labels, id_type vertex, id_type dst_comm) = 0;
    };
    double bgll_vertex_mover_optimizer(const graph& g, vector<id_type>& labels, bgll_objective& book);
    void cda_bgll_generic(const graph&g, const vector<id_type>& init_comms, vector< vector<id_type> >& hier_comms, bgll_objective& book);
    void cda_bgll_modularity(const graph& g, const vector<id_type>& init_comms, vector< vector<id_type> >& hier_comms, double resolution_param);
}

#endif	/* COMMUNITY_H */

