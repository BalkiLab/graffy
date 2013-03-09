/*
 * File:   graph_properties.cpp
 * Author: bharath
 *
 * Created on May 15, 2012, 8:26 AM
 */

#include "graph_properties.h"

void CDLib::get_degree_histogram(const graph& g, vector<id_type>& dist, bool in_degrees) {
    vector<id_type> degrees;
    get_degree_sequence(g, degrees, in_degrees);
    get_discrete_distribution<id_type > (degrees, dist);
}

double CDLib::get_degree_distribution(const graph& g, vector<double>& dist, bool in_degrees) {
    vector<id_type> degree_hist;
    dist.clear();
    double expectation = 0;
    get_degree_histogram(g, degree_hist, in_degrees);
    dist.assign(degree_hist.size(), 0);
    for (id_type i = 0; i < degree_hist.size(); i++) {
        dist[i] = (double) degree_hist[i] / (double) g.get_num_nodes();
        expectation += i * dist[i];
    }
    return expectation;
}

double CDLib::get_excess_degree_distribution(const graph& g, vector<double>& dist, bool in_degrees) {
    dist.clear();
    double expectation = 0;
    double mean = get_degree_distribution(g, dist, in_degrees);
    id_type excess_max = dist.size() - 1;
    for (id_type i = 0; i < excess_max; i++) {
        dist[i] = ((i + 1) * dist[i + 1]) / mean;
        expectation += i * dist[i];
    }
    dist.erase(dist.end() - 1);
    return expectation;
}

double CDLib::get_degree_variance(const graph& g, bool in_degrees) {
    if (g.get_num_nodes() <= 0)
        return 0;
    vector<id_type> degrees;
    get_degree_sequence(g, degrees, in_degrees);
    return variance(degrees);
}

double CDLib::get_degree_assortativity_coefficient(const graph& g, vector<double>& assortativity) {
    /* Gives the assortativity contributions from each node */
    /* The calculation is based on Classifying Complex Networks using Unbiased Local Assortativity, 2010 paper by Piraveenan et al. */
    /* This should only be used for undirected graphs, otherwise this will lead to incorrect results. */
    assortativity.clear();
    double assortativity_coef = 0;
    vector<double> degree_distt;
    double mean = get_degree_distribution(g, degree_distt, 1); // For undirected graph, both in_degrees and out_degrees.
    vector<double> excess_degree_distt(degree_distt.size() - 1, 0);
    double edd_mean = 0, edd_variance = 0;
    for (id_type i = 0; i < excess_degree_distt.size(); i++) {
        excess_degree_distt[i] = ((i + 1) * degree_distt[i + 1]) / mean;
        double imean = i * excess_degree_distt[i];
        edd_mean += imean;
        edd_variance += i*imean;
    }
    edd_variance -= (edd_mean * edd_mean);
    if (edd_variance == 0.0) {
        assortativity.assign(g.get_num_nodes(), ((double) 1 / g.get_num_nodes()));
        return 1;
    } else {
        assortativity.assign(g.get_num_nodes(), 0);
    }
    id_type edge_factor = 2 * g.get_num_edges();
#ifdef ENABLE_MULTITHREADING
        #pragma omp parallel for schedule(dynamic,20) shared(g,assortativity)
#endif
    for(id_type i=0;i<g.get_num_nodes();i++) {
        double avg_excess_degree_neighbour = 0;
        for (adjacent_edges_iterator aeit = g.in_edges_begin(i); aeit != g.in_edges_end(i); aeit++) {
            avg_excess_degree_neighbour += g.get_node_in_degree(aeit->first);
        }
        if (g.get_node_in_degree(i) > 0) {
            avg_excess_degree_neighbour = (avg_excess_degree_neighbour / g.get_node_in_degree(i)) - 1;
            assortativity[i] = ((g.get_node_in_degree(i) - 1) * g.get_node_in_degree(i) * (avg_excess_degree_neighbour - edd_mean)) / (edd_variance * edge_factor);
        }
        assortativity_coef += assortativity[i];
    }
    return assortativity_coef;
}

double CDLib::unbiased_assortativity(const graph& g) {
    /* This is the new way of unbiased assortativity for undirected graph. This is weighted by region-wise assorativity. */
    double weighted_assortativity = 0;
    vector<double> assortativity_node;
    vector<double> degree_distt;
    get_degree_assortativity_coefficient(g, assortativity_node);
    get_degree_distribution(g, degree_distt, 0);
    vector<double> assortativity(degree_distt.size(), 0);
    for (id_type i = 0; i < assortativity_node.size(); i++)
        assortativity[g.get_node_in_degree(i)] += assortativity_node[i];
    for (id_type i = 1; i < degree_distt.size(); i++)
        weighted_assortativity += degree_distt[i] * assortativity[i];
    return weighted_assortativity;
}

double CDLib::difference_assortativity(const graph& g) {
    double harmonic = 0;
#ifdef ENABLE_MULTITHREADING
    #pragma omp parallel for schedule(dynamic,20) shared(g) reduction(+:harmonic)
#endif
    for (id_type i = 0; i < g.get_num_nodes(); i++)
        for (adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i); aeit++)
            harmonic += 1 /(1 +  abs(g.get_node_out_weight(i) - g.get_node_in_weight(aeit->first)));
    return harmonic/(2*g.get_total_weight());
}

double newman_assortativity_directed(const graph& g) {
    double edge_factor = g.get_num_edges();
    double assortaivity_coefficient = 0;
    vector<double> excess_degree_in, excess_degree_out;
    double first_moment_excess_in = 0, second_moment_excess_in = 0, first_moment_excess_out = 0, second_moment_excess_out = 0;
    double mean_degree_in = get_degree_distribution(g, excess_degree_in, 1);
    double mean_degree_out = get_degree_distribution(g, excess_degree_out, 0);
    id_type max_excess_in = excess_degree_in.size() - 1;
    for (id_type i = 0; i < max_excess_in; i++) {
        excess_degree_in[i] = ((i + 1) * excess_degree_in[i + 1]) / mean_degree_in;
        double first = i * excess_degree_in[i];
        first_moment_excess_in += first;
        second_moment_excess_in += i * first;
    }
    excess_degree_in.erase(excess_degree_in.end() - 1);
    double std_excess_in = sqrt(second_moment_excess_in - (first_moment_excess_in * first_moment_excess_in));
    id_type max_excess_out = excess_degree_out.size() - 1;
    for (id_type i = 0; i < max_excess_out; i++) {
        excess_degree_out[i] = ((i + 1) * excess_degree_out[i + 1]) / mean_degree_out;
        double first = i * excess_degree_out[i];
        first_moment_excess_out += first;
        second_moment_excess_out += i * first;
    }
    excess_degree_out.erase(excess_degree_out.end() - 1);
    double std_excess_out = sqrt(second_moment_excess_out - (first_moment_excess_out * first_moment_excess_out));
    vector< vector<double> > joint_distt(excess_degree_out.size(), vector<double>(excess_degree_in.size(), 0));
    for (id_type i = 0; i < g.get_num_nodes(); i++) {
        for (adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i); aeit++) {
            //This place will not be reached when g.get_node_out_degree(i) = 0 OR g.get_node_in_degree(aeit->first) = 0;
            joint_distt[g.get_node_out_degree(i) - 1][g.get_node_in_degree(aeit->first) - 1] += 1 / edge_factor;
        }
    }
#pragma omp parallel for schedule(dynamic,20) shared(g,excess_degree_out,excess_degree_in,joint_distt) reduction(+:assortaivity_coefficient)
    //    Loop from 1 as 0 doesn't add to any contribution to assortativity_coefficient.
    for (id_type j = 1; j < excess_degree_out.size(); j++) {
        for (id_type k = 1; k < excess_degree_in.size(); k++) {
            assortaivity_coefficient += j * k * (joint_distt[j][k]-(excess_degree_out[j] * excess_degree_in[k]));
        }
    }
    assortaivity_coefficient /= (std_excess_out * std_excess_in);
    return assortaivity_coefficient;
}

double newman_assortativity_undirected(const graph& g) {
    /* The calculation is based on Assortative Mixing in Networks, 2002 paper by MEJ Newman */
    double edge_factor = 2 * g.get_num_edges();
    double assortaivity_coefficient = 0;
    vector<double> excess_degree;
    double first_moment_excess = 0, second_moment_excess = 0;
    double mean_degree = get_degree_distribution(g, excess_degree, 0);
    id_type max_excess = excess_degree.size() - 1;
    for (id_type i = 0; i < max_excess; i++) {
        excess_degree[i] = ((i + 1) * excess_degree[i + 1]) / mean_degree;
        double first = i * excess_degree[i];
        first_moment_excess += first;
        second_moment_excess += i * first;
    }
    excess_degree.erase(excess_degree.end() - 1);
    double variance_excess = second_moment_excess - (first_moment_excess * first_moment_excess);
    vector< vector<double> > joint_distt(excess_degree.size(), vector<double>(excess_degree.size(), 0));
    for (id_type i = 0; i < g.get_num_nodes(); i++) {
        for (adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i); aeit++) {
            //This place will not be reached when g.get_node_out_degree(i) = 0 OR g.get_node_out_degree(aeit->first) = 0;
            joint_distt[g.get_node_out_degree(i) - 1][g.get_node_out_degree(aeit->first) - 1] += 1 / edge_factor;
        }
    }
#pragma omp parallel for schedule(dynamic,20) shared(g,excess_degree,joint_distt) reduction(+:assortaivity_coefficient)
    //    Loop from 1 as 0 doesn't add to any contribution to assortativity_coefficient.
    for (id_type j = 1; j < excess_degree.size(); j++) {
        for (id_type k = 1; k < excess_degree.size(); k++) {
            assortaivity_coefficient += j * k * (joint_distt[j][k]-(excess_degree[j] * excess_degree[k]));
        }
    }
    assortaivity_coefficient /= variance_excess;
    return assortaivity_coefficient;
}

double CDLib::get_degree_assortativity_coefficient(const graph& g) {
    /* The calculation is based on Assortative Mixing in Networks, 2002 paper by MEJ Newman */
    if (g.get_num_nodes() <= 0)
        return 0;
    if (g.is_directed())
        return newman_assortativity_directed(g);
    else
        return newman_assortativity_undirected(g);
}

double CDLib::get_rich_club_coefficient(const graph& g, id_type start_hub_degree_def) {
    /* This implementation is accordance to The rich-club phenomenon in the Internet topology, 2004 paper by S. Zhou and R. J. Mondragon */
    id_type num_rich_nodes = 0, num_rich_edges = 0;
#ifdef ENABLE_MULTITHREADING
        #pragma omp parallel for schedule(dynamic,20) shared(g) reduction(+:num_rich_nodes,num_rich_edges)
#endif
    for(id_type i=0;i<g.get_num_nodes();i++) {
        if (g.get_node_out_degree(i) >= start_hub_degree_def) {
            num_rich_nodes++;
            for (adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i); aeit++) {
                if (g.get_node_out_degree(aeit->first) >= start_hub_degree_def) num_rich_edges++;
            }
        }
    }
    return (double) (2 * num_rich_edges) / (num_rich_nodes * (num_rich_nodes - 1));
}

double CDLib::normalized_rich_club_coefficient(const graph& g, id_type start_hub_degree_def) {
    /* This implementation is accordance to Detecting rich-club ordering in complex networks, 2006 paper by Colizza et al. */
    if (start_hub_degree_def <= 0) start_hub_degree_def = 1;
    id_type num_rich_nodes = 0, num_rich_edges = 0;
    double rich_club_uncorrelated = 0, rich_club = 0;
#ifdef ENABLE_MULTITHREADING
        #pragma omp parallel for schedule(dynamic,20) shared(g) reduction(+:num_rich_nodes,num_rich_edges)
#endif
    for(id_type i=0;i<g.get_num_nodes();i++) {
        rich_club_uncorrelated += g.get_node_out_degree(i);
        if (g.get_node_out_degree(i) >= start_hub_degree_def) {
            num_rich_nodes++;
            for (adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i); aeit++) {
                if (g.get_node_out_degree(aeit->first) >= start_hub_degree_def) num_rich_edges++;
            }
        }
    }
    rich_club = (double) (2 * num_rich_edges) / (num_rich_nodes * (num_rich_nodes - 1));
    rich_club_uncorrelated = (start_hub_degree_def * start_hub_degree_def) / rich_club_uncorrelated;
    return rich_club / rich_club_uncorrelated;
}

double CDLib::get_poor_club_coefficient(const graph& g, id_type start_hub_degree_def) {
    /* This implementation is accordance to Detecting rich-club ordering in complex networks, 2006 paper by Colizza et al. */
    id_type num_rich_nodes = 0, num_rich_edges = 0;
#ifdef ENABLE_MULTITHREADING
        #pragma omp parallel for schedule(dynamic,20) shared(g) reduction(+:num_rich_nodes,num_rich_edges)
#endif
    for(id_type i=0;i<g.get_num_nodes();i++) {
        if (g.get_node_out_degree(i) < start_hub_degree_def) {
            num_rich_nodes++;
            for (adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i); aeit++) {
                if (g.get_node_out_degree(aeit->first) < start_hub_degree_def) num_rich_edges++;
            }
        }
    }
    return (double) (2 * num_rich_edges) / (num_rich_nodes * (num_rich_nodes - 1));
}

double CDLib::kl_divergence_from_random_graph(const graph& g) {
    //    Reports symmetric KL-Divergence of the Degree Distribution from Random graph of same size.
    graph er_graph(0, 0);
    generate_erdos_renyi_graph(er_graph, g.get_num_nodes(), g.get_num_edges());
    vector<double> distt1, distt2;
    get_degree_distribution(g, distt1, 0);
    get_degree_distribution(er_graph, distt2, 0);
    return kl_divergence_symmetric(distt1, distt2);
}

double CDLib::distance_from_random_graph(const graph& g, bool hellinger) {
    //    If bool is set, it return Hellinger Distance of the Degree Distribution from Random Graph of same size,
    //    else it returns Bhattacharyya Distance of the Degree Distribution.
    vector<double> distt1, distt2;
    double lambda = (2 * g.get_num_edges()) / g.get_num_nodes();
    id_type degree = 0;
    double probability = exp(-1 * lambda);
    double factorial = 1;
    while (((probability != HUGE_VAL) || (probability != -HUGE_VAL)) && ((degree < lambda) || (probability > 0.001))) {
        distt2.push_back(probability);
        degree++;
        factorial *= degree;
        probability = (pow(lambda, degree) * exp(-1 * lambda)) / factorial;
    }
    get_degree_distribution(g, distt1, 0);
    if (hellinger)
        return hellinger_distance(distt1, distt2);
    else
        return bhattacharyya_distance(distt1, distt2);
}

double CDLib::connectivity_entropy(const graph& g) {
    // Connective Entropy of the Network or Information Entropy of the Network
    double entropy = 0;
    if ((2 * g.get_num_edges()) > 0) {
        for (id_type i = 0; i < g.get_num_nodes(); i++) {
            double prob = (double) g.get_node_out_degree(i) / (2 * g.get_num_edges());
            if (prob != 0) {
                entropy += prob * (log(prob) / log(2));
            }
        }
        entropy *= -1;
    }
    entropy /= log(g.get_num_nodes()) / log(2); // Normalization of Entropy
    return entropy;
}

double CDLib::graph_modularity(graph& g) {
    /* Functions returns the modularity of the whole graph considering the whole graph as a single community */
    double denominator = 2 * g.get_num_edges();
    double sum = 0;
    if (g.is_directed()) {
        for (id_type i = 0; i < g.get_num_nodes(); i++) {
            for (id_type j = 0; j < g.get_num_nodes(); j++) {
                sum += g.get_edge_weight(i, j) - ((g.get_node_out_degree(i) * g.get_node_out_degree(j)) / denominator);
            }
        }
    } else {
        for (id_type i = 0; i < g.get_num_nodes(); i++) {
            for (id_type j = i + 1; j < g.get_num_nodes(); j++) {
                sum += 2 * (g.get_edge_weight(i, j) - ((g.get_node_out_degree(i) * g.get_node_out_degree(j)) / denominator));
            }
        }
        for (id_type i = 0; i < g.get_num_nodes(); i++) {
            sum += g.get_edge_weight(i, i) - ((g.get_node_out_degree(i) * g.get_node_out_degree(i)) / denominator);
        }
    }
    return sum;
}
