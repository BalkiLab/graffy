/*
 * File:   epidemic.cpp
 * Author: sudip
 *
 * Created on Jan 4, 2015, 8:00 PM
 */

#include "epidemic.h"

using namespace CDLib;

void CDLib::diffusion(const graph& g, double transition_prob, id_type monte_carlo, vector<double>& cover) {
    cover.clear();
    if (transition_prob < 0)
        return;
    if (monte_carlo < 1)
        monte_carlo = 1;
    srand(time(NULL));
    RandomGenerator<id_type> p_gen(0, g.get_num_nodes() - 1, 0);
    for (id_type i = 0; i < monte_carlo; i++) {
        vector<id_type> step;
        propagation(g, p_gen.next(), transition_prob, step);
        if (step.size() > cover.size()) {
            for (id_type j = cover.size(); j < step.size(); j++)
                cover.push_back(0);
        }
        for (id_type j = 0; j < step.size(); j++)
            cover[j] += step[j];
    }
    for (id_type i = 0; i < cover.size(); i++)
        cover[i] /= monte_carlo;
}

void CDLib::propagation(const graph& g, id_type seed_node_id, double transition_prob, vector<id_type>& step) {
    step.clear();
    if ((seed_node_id >= g.get_num_nodes()) || (transition_prob < 0)) 
        return;
    RandomGenerator<double> p_gen(0, 1, 0);
    vector< vector<id_type> > nodes;
    node_set visited;
    id_type counts = 1;
    nodes.push_back(vector<id_type>());
    nodes[0].push_back(seed_node_id);
    visited.insert(seed_node_id);
    while (nodes.size() <= counts) {
        nodes.push_back(vector<id_type>());
        for (id_type i = 0; i < nodes[counts - 1].size(); i++) {            
            for (adjacent_edges_iterator aeit = g.out_edges_begin(nodes[counts - 1][i]); aeit != g.out_edges_end(nodes[counts - 1][i]); aeit++) {
                if ((p_gen.next() <= transition_prob) && (visited.find(aeit->first) == visited.end())) {
                    visited.insert(aeit->first);
                    nodes[counts].push_back(aeit->first);
                }
            }
        }
        if (nodes[counts].size() != 0)
            counts++;
    }
    for (id_type i = 0; i < nodes.size(); i++) {
        step.push_back(nodes[i].size());
    }
    step.pop_back();
}



