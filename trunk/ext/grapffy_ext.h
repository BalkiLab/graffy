/*
 * File:   grapffy_ext.h
 * Author: sudip
 *
 * Created on 2 March, 2013, 11:11 PM
 */

#ifndef GRAPFFY_EXT_H
#define	GRAPFFY_EXT_H

#include "includes.h"
#include "R_ext.h"
#include "igraph_ext.h"
#include "octave_ext.h"

struct graph_assortativity {
    vector<double> assortativity_distt, degree_assortativity;
    double assortaivity_coefficient, head_assortativity, tail_assortativity, max_assortative_cov;
    double rich_club_coefficient_from_xmin, poor_club_below_xmin;

    graph_assortativity() {
        assortaivity_coefficient = 0;
        head_assortativity = 0;
        tail_assortativity = 0;
        max_assortative_cov = 0;
        rich_club_coefficient_from_xmin = 0;
        poor_club_below_xmin = 0;
        assortativity_distt.clear();
        degree_assortativity.clear();
    }

    void local_assortativity(const graph& g, const vector<double>& degree_distt, const statistics<id_type>& deg_stat, const struct fitresult_R & degree_fit) {
        /* The calculation is based on Classifying Complex Networks using Unbiased Local Assortativity, 2010 paper
         * by Piraveenan et al. and also rich club coefficient based on Detecting rich-club ordering in complex
         * networks, 2006 paper by Colizza et al.*/
        vector<double> excess_degree_distt(deg_stat.max_val, 0);
        double edd_mean = 0, edd_variance = 0;
        for (id_type i = 0; i < excess_degree_distt.size(); i++) {
            excess_degree_distt[i] = ((i + 1) * degree_distt[i + 1]) / deg_stat.mean_val;
            double imean = i * excess_degree_distt[i];
            edd_mean += imean;
            edd_variance += i*imean;
        }
        edd_variance -= (edd_mean * edd_mean);
        if (edd_variance == 0) {
            assortativity_distt.assign(g.get_num_nodes(), (1 / g.get_num_nodes()));
            assortaivity_coefficient = 1;
            return;
        } else {
            assortativity_distt.assign(g.get_num_nodes(), 0);
        }
        vector<double>& loc_assortativity_distt = assortativity_distt;
        double loc_assortaivity_coefficient = 0;
        id_type num_rich_nodes = 0, num_rich_edges = 0, num_poor_edges = 0, hub_deg = degree_fit.xmin;
        id_type edge_factor = 2 * g.get_num_edges();
#pragma omp parallel for schedule(dynamic,20) shared(g,loc_assortativity_distt) reduction(+:loc_assortaivity_coefficient)
        for (id_type i = 0; i < g.get_num_nodes(); i++) {
            double avg_excess_degree_neighbour = 0;
            if (g.get_node_out_degree(i) >= hub_deg) num_rich_nodes++;
            for (adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i); aeit++) {
                avg_excess_degree_neighbour += g.get_node_out_degree(aeit->first);
                if ((g.get_node_out_degree(i) >= hub_deg) && (g.get_node_out_degree(aeit->first) >= hub_deg)) num_rich_edges++;
                if ((g.get_node_out_degree(i) < hub_deg) && (g.get_node_out_degree(aeit->first) < hub_deg)) num_poor_edges++;
            }
            if (g.get_node_out_degree(i) > 0) {
                avg_excess_degree_neighbour = (avg_excess_degree_neighbour / g.get_node_out_degree(i)) - 1;
                loc_assortativity_distt[i] = ((g.get_node_out_degree(i) - 1) * g.get_node_out_degree(i) * (avg_excess_degree_neighbour - edd_mean)) / (edd_variance * edge_factor);
                loc_assortaivity_coefficient += loc_assortativity_distt[i];
            }
        }
        rich_club_coefficient_from_xmin = (double) num_rich_edges / (num_rich_nodes * (num_rich_nodes - 1));
        poor_club_below_xmin = (double) num_poor_edges / ((g.get_num_nodes() - num_rich_nodes) * ((g.get_num_nodes() - num_rich_nodes) - 1));
        assortaivity_coefficient = loc_assortaivity_coefficient;
    }

    void regional_assortativity(const graph& g, const statistics<id_type>& deg_stat, const struct fitresult_R & degree_fit) {
        degree_assortativity.assign(deg_stat.max_val + 1, 0);
        for (unsigned long i = 0; i < assortativity_distt.size(); i++) {
            degree_assortativity[g.get_node_out_degree(i)] += assortativity_distt[i];
            if (degree_fit.gof_tail_pl < numeric_limits<double>::infinity()) {
                if (g.get_node_out_degree(i) < degree_fit.xmin) {
                    head_assortativity += assortativity_distt[i];
                } else {
                    tail_assortativity += assortativity_distt[i];
                }
            }
        }
    }
};

struct graph_measures {
    const graph& g;
    id_type number_of_nodes, number_of_edges;
    vector<id_type> degrees, dhist;
    statistics<id_type> deg_stat;
    vector<double>degree_distt;
    double density, degree_entropy;

    double diameter, path_entropy, efficiency, fraction_in_lcc, reachibility, avg_path_length_lcc;
    double avg_clustering_coefficient;
    vector<double> distance_distt;

    double distance_from_er;
    struct fitresult_R degree_fit;
    double plfit_frac_in_tail;
    struct graph_assortativity assortativity;
    double dist_degree_distt, dist_distance_distt;

    graph_measures(const graph & org) : g(org) {
        number_of_nodes = 0;
        number_of_edges = 0;
        density = 0;
        diameter = 0;
        degree_entropy = 0;
        path_entropy = 0;
        efficiency = 0;
        distance_from_er = 0;
        reachibility = 0;
        distance_from_er = 0;
        fraction_in_lcc = 0;
        avg_path_length_lcc = 0;
        avg_clustering_coefficient = 0;
        plfit_frac_in_tail = 0;
        dist_degree_distt = 0;
        dist_distance_distt = 0;
    }

    void eval_all_measure() {
        /* This sequence of execution should be maintained. Only this method should be executed from outside. */
        if (g.get_num_nodes() < 1)
            return;
        if (g.get_num_edges() == 0) {
            diameter = numeric_limits<double>::infinity();
            fraction_in_lcc = (double) 1 / g.get_num_nodes();
            dist_degree_distt = 1;
            dist_distance_distt = 1;
            distance_from_er = 1;
            return;
        }
        eval_deg_seq();
        eval_deg_distt();
        eval_paths(g, distance_distt);
        eval_er_distance();
        eval_plfit();
        eval_assortaivity();
    }

    void eval_all_measure(graph_measures & gm) {
        /* This sequence of execution should be maintained. Only this method should be executed from outside. */
        if (g.get_num_nodes() == 0)
            return;
        eval_all_measure();
        eval_compare_distribution(gm);
    }

private:

    void eval_deg_seq() {
        number_of_nodes = g.get_num_nodes();
        number_of_edges = g.get_num_edges();
        density = g.get_density();
        get_degree_sequence(g, degrees, 0); // 0: works on out_degrees of nodes
        compute_statistics<id_type > (degrees, deg_stat);
    }

    void eval_deg_distt() {
        dhist.assign(deg_stat.max_val + 1, 0);
        for (id_type i = 0; i < degrees.size(); i++)
            dhist[degrees[i]]++;
        degree_distt.assign(dhist.size(), 0);
        for (id_type i = 0; i < dhist.size(); i++) {
            degree_distt[i] = (double) dhist[i] / (double) g.get_num_nodes();
            if (degree_distt[i] != 0)
                degree_entropy -= degree_distt[i] * (log(degree_distt[i]) / log(2));
        }
    }

    void eval_er_distance() {
        graph er_graph(0, 0);
        generate_erdos_renyi_graph(er_graph, g.get_num_nodes(), g.get_num_edges());
        vector<double> er_distt;
        get_degree_distribution(er_graph, er_distt, 0);
        distance_from_er = hellinger_distance(degree_distt, er_distt);
    }

    void eval_paths(const graph& loc_g, vector<double>& loc_distance_distt) {
        /* The code works properly only for undirected graph. */
        vector<double> reachability_distt(g.get_num_nodes(), 0);
        vector<double> total_distance(g.get_num_nodes(), 0.0);
        loc_distance_distt.assign(g.get_num_nodes() + 1, 0.0);
        id_type all_sources = 0;
        double loc_diameter = 0;
        bool is_diameter_infinite = false;
        double loc_eff = 0, loc_cc = 0;
        double normalization = g.get_num_nodes() * (g.get_num_nodes() - 1);
        vector<bool> has_node_visited(g.get_num_nodes(), 0);
        vector<node_set> components;
#pragma omp parallel for schedule(dynamic,20) shared(loc_g,loc_distance_distt,loc_diameter,total_distance,reachability_distt,has_node_visited,components) reduction(+:all_sources,loc_eff,loc_cc)
        for (id_type i = 0; i < loc_g.get_num_nodes(); i++) {
            node_set current_comm;
            current_comm.insert(i);
            unordered_map<id_type, id_type> local_dist;
            vector<double> distances(loc_g.get_num_nodes(), numeric_limits<double>::infinity());
            queue<id_type> q_bfs;
            q_bfs.push(i);
            distances[i] = 0;
            while (!q_bfs.empty()) {
                id_type current = q_bfs.front();
                q_bfs.pop();
                for (adjacent_edges_iterator aeit = loc_g.out_edges_begin(current); aeit != loc_g.out_edges_end(current); aeit++) {
                    if (distances[aeit->first] == numeric_limits<double>::infinity() && current != aeit->first) {
                        distances[aeit->first] = distances[current] + 1;
                        q_bfs.push(aeit->first);
                        pair < unordered_map<id_type, id_type>::iterator, bool> ret = local_dist.insert(make_pair(distances[aeit->first], 1));
                        if (!ret.second) ret.first->second++;
                        loc_eff += (1 / distances[aeit->first]);
                        total_distance[i] += distances[aeit->first];
                        reachability_distt[i]++;
                        all_sources++;
                        current_comm.insert(aeit->first);
                    }
                }
            }
            loc_cc += node_clustering_coefficient(loc_g, i);
            local_dist.insert(make_pair(loc_g.get_num_nodes(), loc_g.get_num_nodes() - reachability_distt[i] - 1));
#pragma omp critical
            {
                for (unordered_map<id_type, id_type>::iterator it = local_dist.begin(); it != local_dist.end(); it++) {
                    loc_distance_distt[it->first] += it->second;
                    if (it->first != loc_g.get_num_nodes() && it->first > loc_diameter) loc_diameter = it->first;
                    if (it->first == loc_g.get_num_nodes() && it->second) is_diameter_infinite = true;
                }
                if (!has_node_visited[i]) {
                    components.push_back(current_comm);
                    has_node_visited[i] = true;
                    for (node_set::iterator it = current_comm.begin(); it != current_comm.end(); it++) {
                        has_node_visited[*it] = true;
                    }
                }
            }

        }
        distance_distt[loc_diameter + 4] = distance_distt[g.get_num_nodes()]; // Allowing 3 zero values for better plotting
        distance_distt.erase(distance_distt.begin() + loc_diameter + 5, distance_distt.end());
        if (is_diameter_infinite) {
            diameter = numeric_limits<double>::infinity();
        } else {
            diameter = loc_diameter;
        }
        loc_eff /= normalization;
        avg_clustering_coefficient = loc_cc / g.get_num_nodes();
        efficiency = loc_eff;
        /* This consider only when there is atleast one edge. */
        if (all_sources > 0) {
            vector<id_type> index_max;
            id_type lcc = 0;
            for (id_type i = 0; i < g.get_num_nodes(); i++) {
                reachibility += reachability_distt[i];
                double reach_entropy = reachability_distt[i] / all_sources;
                if (reach_entropy != 0) {
                    path_entropy -= reach_entropy * (log(reach_entropy) / log(2));
                }
                if (i < components.size()) {
                    if (components[i].size() >= lcc) {
                        if (components[i].size() > lcc) {
                            index_max.clear();
                            lcc = components[i].size();
                            index_max.push_back(i);
                        } else {
                            index_max.push_back(i);
                        }
                    }
                }
            }
            reachibility /= normalization;
            fraction_in_lcc = (double) lcc / g.get_num_nodes();
            /* Any measures for LCC can be calculated in the below loop */
            for (id_type i = 0; i < index_max.size(); i++) {
                for (node_set::iterator it = components[index_max[i]].begin(); it != components[index_max[i]].end(); it++) {
                    avg_path_length_lcc += total_distance[*it];
                }
            }
            if (index_max.size() > 0) {
                avg_path_length_lcc /= (lcc * (lcc - 1) * index_max.size());
            }
        }
        for (unsigned long i = 0; i < distance_distt.size(); i++) {
            distance_distt[i] /= normalization;
        }
    }

    void eval_plfit() {
        degree_fit.fit_degree_distribution(g, 1); // 0: works on out_degrees of nodes.
        for (id_type i = 0; i < degree_distt.size(); i++) {
            if (i >= degree_fit.xmin) {
                plfit_frac_in_tail += degree_distt[i];
            }
        }
    }

    void eval_assortaivity() {
        assortativity.local_assortativity(g, degree_distt, deg_stat, degree_fit);
        assortativity.regional_assortativity(g, deg_stat, degree_fit);
    }

    void eval_compare_distribution(graph_measures & gm) {
        dist_degree_distt = hellinger_distance(degree_distt, gm.degree_distt);
        dist_distance_distt = hellinger_distance(distance_distt, gm.distance_distt);
    }

    void display_stats() {
        cout << "     Displaying STATS of " << g.get_graph_name() << endl;
        cout << "--------------------------------------------------" << endl;
        cout << "Nodes: " << number_of_nodes << "\t\tEdges: " << number_of_edges << "\t\tDensity: " << density << endl;
        cout << "Degree Stats :" << endl;
        cout << "Min Degree: " << deg_stat.min_val << "\t\tMax Degree: " << deg_stat.max_val << endl;
        cout << "Avg. Degree: " << deg_stat.mean_val << "\t\tVariance: " << deg_stat.variance << endl;
    }
};

struct graph_spectral_properties {
    double spectral_radius_adjacency, spectral_gap_adjacency, largest_evalue_laplacian;

    graph_spectral_properties() {
        spectral_radius_adjacency = 0; // Largest EigenValue of Graph Adjacency Matrix
        spectral_gap_adjacency = 0; // Difference between Largest and Second Largest EigenValue of Graph Adjacency Matrix
        largest_evalue_laplacian = 0; // Largest EigenValue of Graph Laplacian Matrix
    }

    graph_spectral_properties(graph & g) {
        graph_spectral_properties();
        get_spectral_properties(g);
    }

    void get_spectral_properties(graph & g) {
        string result = get_spectral(g);
        stringstream ss(result);
        string name;
        ss >> name >> spectral_radius_adjacency >> spectral_gap_adjacency >> largest_evalue_laplacian;
    }
};

double modularity_blondel(const graph& g) {
    igraph_ext *gcom = new igraph_ext(g);
    gcom->blondel_multilevel();
    vector<node_set> comms;
    gcom->return_community_structure(comms);
    delete gcom;
    return modularity(g, comms);
}

#endif	/* GRAPFFY_EXT_H */

