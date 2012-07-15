/* 
 * File:   matrix.h
 * Author: sudip
 *
 * Created on 22 May, 2012, 10:21 PM
 */

#ifndef MATRIX_H
#define	MATRIX_H



#define EIGEN_USE_MKL_ALL
#include <eigen3/Eigen/Eigen>
using namespace Eigen;
#include "graph.h"

namespace CDLib {
    void eigen_decomposition_of_laplacian(graph& g,vector<double>& evals, vector< vector<double> >& evecs);
    void eigen_decomposition_of_adjacency(graph& g, vector< std::complex<double> >& evals, vector< vector< std::complex<double> > >& evecs);
    void resistance_distance(graph& g, vector< vector<double> >& resistance);
    void conductance_distance(graph& g, vector< vector<double> >& conductance);
    void hop_transition_probability(graph& g, vector< vector<double> >& hop_transition_probability);
    void average_distance(graph& g, vector< vector<double> >& hop_distance);
};

#endif	/* MATRIX_H */

