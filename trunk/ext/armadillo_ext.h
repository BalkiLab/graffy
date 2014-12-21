/*
 * File:   armadillo_ext.h
 * Author: bharath
 *
 * Created on 1 March, 2013, 10:18 PM
 */

#ifndef ARMADILLO_EXT_H
#define	ARMADILLO_EXT_H

#include "includes.h"
#include <armadillo>
using namespace std;
namespace armadillo {

    template <typename T>
    void adjacency_matrix(const CDLib::graph& g, Mat<T>& A) {
        size_t mat_sz = g.get_num_nodes();
        A.zeros(mat_sz, mat_sz);
        for (CDLib::id_type i = 0; i < mat_sz; i++)
            for (CDLib::adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i); aeit++)
                if (aeit->first != i)
                    A(i, aeit->first) = aeit->second;
    }

    template <typename T>
    void degree_matrix(const CDLib::graph& g, Mat<T>& D) {
        size_t mat_sz = g.get_num_nodes();
        D.zeros(mat_sz, mat_sz);
        for (CDLib::id_type i = 0; i < mat_sz; i++)
            D(i, i) = g.get_node_out_weight(i);
    }

    template <typename T>
    void laplacian_matrix(const CDLib::graph& g, Mat<T>& L) {
        size_t mat_sz = g.get_num_nodes();
        L.zeros(mat_sz, mat_sz);
        for (CDLib::id_type i = 0; i < mat_sz; i++) {
            for (CDLib::adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i); aeit++) {
                if (aeit->first != i) L(i, aeit->first) = -aeit->second;
                else L(i, i) = g.get_node_out_weight(i);

            }
        }
    }

    template <typename T>
    void transition_probability_matrix(const CDLib::graph& g, Mat<T>& P) {
        size_t mat_sz = g.get_num_nodes();
        P.zeros(mat_sz, mat_sz);
        for (CDLib::id_type i = 0; i < mat_sz; i++)
            for (CDLib::adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i); aeit++)
                if (aeit->first != i)
                    P(i, aeit->first) = aeit->second / g.get_node_out_weight(i);
    }

    template <typename T>
    void transition_probability_matrix_botgrep(const CDLib::graph& g, Mat<T>& P) {
        size_t mat_sz = g.get_num_nodes();
        P.zeros(mat_sz, mat_sz);
        for (CDLib::id_type i = 0; i < mat_sz; i++)
            for (CDLib::adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i); aeit++)
                if (aeit->first != i)
                    P(i, aeit->first) = min(aeit->second / g.get_node_out_weight(i), aeit->second / g.get_node_out_weight(aeit->first));
    }

}


#endif	/* ARMADILLO_EXT_H */

