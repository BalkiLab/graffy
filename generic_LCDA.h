/* 
 * File:   generic_LCDA.h
 * Author: prashant
 *
 * Created on July 11, 2012, 5:51 PM
 */

#ifndef GENERIC_LCDA_H
#define	GENERIC_LCDA_H

#include "graph.h"
using namespace std;

namespace CDLib {

    struct lcd_core {
        node_set C;
        node_set B;
        node_set U;
        double I;
        double E;

        lcd_core() : C(), B(), U(), I(0), E(0) {
        }

        lcd_core(const graph& g, id_type seed) : C(), B(), U(), I(0), E(0) {
            C.insert(seed);
            B.insert(seed);
            for(adjacent_edges_iterator aeit = g.out_edges_begin(seed); aeit != g.out_edges_end(seed); aeit++)
                U.insert(aeit->first);
            E = g.get_node_out_degree(seed);
            //I=0; ??????????????????????????
        }
    };

    struct lcd_new {
        double I_new;
        double E_new;
        double n_new;
        double n_b_new;
        double n_u_new;

        lcd_new() : I_new(0), E_new(0), n_new(0), n_b_new(0), n_u_new(0) {
        }
    };

    typedef double (*lcd_obj_fn)(const graph&, lcd_new&, id_type);
    typedef bool (*lcd_stop_fn)(const graph&, lcd_core&, node_set&, node_set&, vector<double>&); //g and deleteQ not needed till now
    typedef void (*lcd_traversal_fn)(const graph&, lcd_core&, lcd_obj_fn, node_set&);

    void local_community_detection_std(
            const graph& g,
            id_type seed,
            vector<double>& stop_fn_params,
            lcd_traversal_fn traversal,
            bool include_deletion_phase,
            lcd_obj_fn,
            lcd_obj_fn,
            lcd_stop_fn,
            node_set& comm);

};

#endif	/* GENERIC_LCDA_H */
