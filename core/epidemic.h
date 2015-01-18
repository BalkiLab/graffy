/* 
 * File:   epidemic.h
 * Author: sudip
 *
 * Created on 4 January, 2015, 8:00 PM
 */

#ifndef EPIDEMIC_H
#define	EPIDEMIC_H

#include "graph.h"
#include "datastructures.h"
using namespace std;

namespace CDLib {
    void propagation(const graph& g, id_type seed_node_id, double transition_prob, vector<id_type>& step);
    void diffusion_cover(const graph& g, double transition_prob, id_type monte_carlo, vector<double>& cover);
    void diffusion_step(const graph& g, double transition_prob, id_type monte_carlo, vector<double>& cover);
};

#endif	/* EPIDEMIC_H */

