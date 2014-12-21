/*
 * File:   eigen_ext.h
 * Author: bharath
 *
 * Created on 1 March, 2013, 10:19 PM
 */

#ifndef EIGEN_EXT_H
#define	EIGEN_EXT_H

#define EIGEN_USE_MKL_ALL
#include "includes.h"

#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Dense>

struct pca_with_eigen {
    /* evals and evecs are the eigenvalue and eigenvector of the covariance matrix. */
    vector<double> evals, leading_component, cummulative_variance;
    vector < vector<double > evecs, covariance_matrix, all_components, components;

    pca_with_eigen(vector< vector<double> >& features) {
        /* Each row of the matrix contains a feature. Each column is a data point. */
        /* The original feature value will get changed and get adjusted to zero mean. */
        adjust_mean_feature_matrix(features);
        make_covariance_matrix(features);
        eigen_decomposition_of_matrix();
        get_increasing_variance();
        principal_component(features);
    }

    /* This will give components upto the give position. Will also include the leading component in 0th position. */
    void populate_components(vector< vector<double> >& features, unsigned long number_of_components) {
        if ((evals.size() == 0) || (features.size() == 0))
            return;
        if (number_of_components > evals.size()) {
            perror("Number of components exceed number of features.");
            return;
        }
        components.clear();
        components.assign(number_of_components, vector<double>());
        components[0].assign(leading_component.begin(), leading_component.end());
        unsigned long change = evals.size() - 1; // Will help to arrange in decreasing order from increasing order.
        for (unsigned long j = 0; j < features[0].size(); j++) {
            for (id_type x = 1; x < number_of_components; x++) {
                components[j].assign(evals.size(), 0);
                for (id_type k = 0; k < evals.size(); k++) {
                    components[j][x] += features[k][j] * evecs[k][change - x];
                }
            }
        }
    }

    /* This function needs to be executed separately to get all components after transformation */
    void populate_all_components(vector< vector<double> >& features) {
        if ((evals.size() == 0) || (features.size() == 0))
            return;
        components.clear();
        components.assign(evals.size(), vector<double>());
        components[0].assign(leading_component.begin(), leading_component.end());
        unsigned long change = evals.size() - 1; // Will help to arrange in decreasing order from increasing order.
        for (unsigned long j = 0; j < features[0].size(); j++) {
            for (id_type x = 1; x < evals.size(); x++) {
                all_components[j].assign(evals.size(), 0);
                for (id_type k = 0; k < evals.size(); k++) {
                    all_components[j][x] += features[k][j] * evecs[k][change - x];
                }
            }
        }
    }
private:

    void adjust_mean(vector<double>& data) {
        if (data.size() < 1)
            return;
        double mean = accumulate(data.begin(), data.end(), 0.0);
        mean /= data.size();
        for (id_type i = 0; i < data.size(); i++) {
            data[i] -= mean;
        }
    }

    void adjust_mean_feature_matrix(vector< vector<double> >& features) {
        if (features.size()) {
            for (id_type i = 0; i < features.size(); i++)
                adjust_mean(features[i]); // Adjust mean of each features represented by rows.
        }
    }

    void make_covariance_matrix(const vector< vector<double> >& features) {
        covariance_matrix.clear();
        if (features.size() == 0)
            return;
        if (features.size() == 1) {
            covariance_matrix.assign(1, vector<double>(1, 0));
            return;
        }
        covariance_matrix.assign(features.size(), vector<double>(features.size(), 0));
        for (id_type i = 0; i < features.size(); i++) {
            for (id_type j = i + 1; j < features.size(); j++) {
                covariance_matrix[i][j] = covariance(features[i], features[j]);
                covariance_matrix[j][i] = covariance_matrix[i][j];
            }
        }
    }

    void eigen_decomposition_of_matrix() {
        // Populates the Eigen Decomposition of the Covariance Matrix given in the form of STL 2D-Vector
        MatrixXd mat;
        MatrixXd_from_double_vector(covariance_matrix, mat);
        VectorXd lambda;
        MatrixXd V;
        eigen_decomposition_SelfAdjoint(mat, lambda, V);
        evals.clear();
        evecs.clear();
        double_vector_from_MatrixXd(V, evecs); // Each column represents each Eigen Vector.
        double_vector_from_VectorXd(lambda, evals);
    }

    void get_increasing_variance() {
        double sum = accumulate(evals.begin(), evals.end(), 0);
        cummulative_variance.assign(evals.size(), 0);
        cummulative_variance[0] = evals[evals.size() - 1] / sum;
        for (unsigned long i = 1, k = evals.size() - 2; k >= 0; k--, i++) {
            cummulative_variance[i] = (cummulative_variance[i - 1] + evals[k]) / sum;
        }
    }

    void principal_component(vector< vector<double> >& features) {
        if ((evals.size() == 0) || (features.size() == 0))
            return;
        int max_id = evals.size() - 1;
        for (unsigned long j = 0; j < features[0].size(); j++) {
            double value = 0;
            for (unsigned long k = 0; k < evals.size(); k++) {
                value += features[k][j] * evecs[k][max_id];
            }
            if (evecs[0][max_id] < 0)
                value *= -1; // Changed the sign to make it consistent across different runs. It will always start with positive value.
            leading_component.push_back(value);
        }
    }
};

void make_eigen_format_adjacency_matrix(graph& g, MatrixXd& gmat) {
    // Creates the adjacency matrix for the given graph
    gmat.resize(g.get_num_nodes(), g.get_num_nodes());
    if (g.is_directed()) {
        for (id_type i = 0; i < g.get_num_nodes(); i++)
            for (id_type j = 0; j < g.get_num_nodes(); j++)
                gmat(i, j) = g.get_edge_weight(i, j);
    } else {
        for (id_type i = 0; i < g.get_num_nodes(); i++) {
            for (id_type j = 0; j < g.get_num_nodes(); j++) {
                gmat(i, j) = g.get_edge_weight(i, j);
                gmat(j, i) = g.get_edge_weight(j, i);
            }
        }
    }
}

void make_eigen_format_laplacian_matrix(graph& g, MatrixXd& lmat) {
    // Creates the laplacian matrix for the given graph
    lmat.resize(g.get_num_nodes(), g.get_num_nodes());
    if (g.is_directed()) {
        for (id_type i = 0; i < g.get_num_nodes(); i++) {
            for (id_type j = 0; j < g.get_num_nodes(); j++) {
                if (g.get_edge_weight(i, j))
                    lmat(i, j) = -1;
                else
                    lmat(i, j) = 0;
            }
        }
    } else {
        for (id_type i = 0; i < g.get_num_nodes(); i++) {
            for (id_type j = 0; j < g.get_num_nodes(); j++) {
                if (g.get_edge_weight(i, j)) {
                    lmat(i, j) = -1;
                    lmat(j, i) = -1;
                } else {
                    lmat(i, j) = 0;
                    lmat(j, i) = 0;
                }
            }
        }
    }
    for (id_type i = 0; i < g.get_num_nodes(); i++)
        lmat(i, i) = g.get_node_out_degree(i);
}

void eigen_decomposition(MatrixXd& mat, VectorXcd& lambda, MatrixXcd& V) {
    // Creates the eigen decomposition of the given matrix into eigen values and eigen vectors
    if (mat.rows() == mat.cols()) {
        EigenSolver<MatrixXd> es(mat);
        lambda = es.eigenvalues();
        V = es.eigenvectors();
    }
}

void eigen_decomposition_SelfAdjoint(MatrixXd& mat, VectorXd& lambda, MatrixXd& V) {
    // Creates the eigen decomposition of the given self adjoint matrix into eigen values and eigen vectors
    if (mat.rows() == mat.cols()) {
        SelfAdjointEigenSolver<MatrixXd> es(mat);
        lambda = es.eigenvalues();
        V = es.eigenvectors();
    }
}

void complex_vector_from_VectorXcd(VectorXcd& v, vector< std::complex<double> >& vecs) {
    // Convert VectorXcd type in Eigen Library to Vector of Complex Double type of C++ STL
    long size = v.size();
    for (long i = 0; i < size; i++)
        vecs.push_back(v(i));
}

void complex_vector_from_MatrixXcd(MatrixXcd& v, vector< vector< std::complex<double> > >& mats) {
    // Convert MatrixXcd type in Eigen Library to Vector of Vector of Complex  Double type of C++ STL
    // Each Column represents 1 Eigen vector
    long rows = v.rows();
    for (long i = 0; i < rows; i++) {
        vector< std::complex<double> > temp;
        VectorXcd operation_vs = v.row(i);
        complex_vector_from_VectorXcd(operation_vs, temp);
        mats.push_back(temp);
    }
}

void double_vector_from_VectorXcd(VectorXcd& v, vector<double>& vecs) {
    // Convert VectorXcd type in Eigen Library to Vector of Double type of C++ STL
    long size = v.size();
    for (long i = 0; i < size; i++)
        vecs.push_back(real(v(i)));
}

void double_vector_from_MatrixXcd(MatrixXcd& v, vector< vector<double> >& mats) {
    // Convert MatrixXcd type in Eigen Library to Vector of Vector of Double type of C++ STL
    // Each Column represents 1 Eigen vector
    long rows = v.rows();
    for (long i = 0; i < rows; i++) {
        vector<double> temp;
        VectorXcd operation_vs = v.row(i);
        double_vector_from_VectorXcd(operation_vs, temp);
        mats.push_back(temp);
    }
}

void double_vector_from_VectorXd(VectorXd& v, vector<double>& vecs) {
    // Convert VectorXd type in Eigen Library to Vector of Double type of C++ STL
    long size = v.size();
    for (long i = 0; i < size; i++)
        vecs.push_back(v(i));
}

void double_vector_from_MatrixXd(MatrixXd& v, vector< vector<double> >& mats) {
    // Convert MatrixXcd type in Eigen Library to Vector of Vector of Double type of C++ STL
    // Each Column represents 1 Eigen vector
    long rows = v.rows();
    for (long i = 0; i < rows; i++) {
        vector<double> temp;
        VectorXd operation_vs = v.row(i);
        double_vector_from_VectorXd(operation_vs, temp);
        mats.push_back(temp);
    }
}

void MatrixXd_from_double_vector(vector< vector<double> >& mats, MatrixXd &m) {
    if ((mats.size() > 0) && (mats.size() == mats[0].size())) {
        m.resize(mats.size(), mats.size());
        for (id_type i = 0; i < mats.size(); i++)
            for (id_type j = 0; j < mats[i].size(); j++)
                m(i, j) = mats[i][j];
    }
}

void VectorXd_from_double_vector(vector<double>& vecs, VectorXd& v) {
    for (id_type i = 0; i < vecs.size(); i++)
        v(i) = vecs[i];
}

void eigen_decomposition_of_matrix(vector< vector<double> >& matrix, vector<double>& evals, vector< vector<double> >& evecs) {
    // Returns the Eigen Decomposition of the Matrix given in the form of STL 2D-Vector
    MatrixXd mat;
    MatrixXd_from_double_vector(matrix, mat);
    VectorXcd lambda;
    MatrixXcd V;
    eigen_decomposition(mat, lambda, V);
    evals.clear();
    evecs.clear();
    /* Eigenvalues are always arranged in ascending order. */
    double_vector_from_MatrixXcd(V, evecs); // Each column represents each Eigen Vector to that Eigenvalue.
    double_vector_from_VectorXcd(lambda, evals);
}

void eigen_decomposition_of_sefadjoint_matrix(vector< vector<double> >& matrix, vector<double>& evals, vector< vector<double> >& evecs) {
    // Returns the Eigen Decomposition of the Matrix given in the form of STL 2D-Vector
    MatrixXd mat;
    MatrixXd_from_double_vector(matrix, mat);
    VectorXd lambda;
    MatrixXd V;
    eigen_decomposition(mat, lambda, V);
    evals.clear();
    evecs.clear();
    /* Eigenvalues are always arranged in ascending order. */
    double_vector_from_MatrixXd(V, evecs); // Each column represents each Eigen Vector to that Eigenvalue.
    double_vector_from_VectorXd(lambda, evals);
}

void eigen_decomposition_of_laplacian(graph& g, vector<double>& evals, vector< vector<double> >& evecs) {
    // Returns the Eigen Decomposition of the Graph Laplacian
    MatrixXd lmat;
    make_eigen_format_laplacian_matrix(g, lmat);
    VectorXcd lambda;
    MatrixXcd V;
    eigen_decomposition(lmat, lambda, V);
    evals.clear();
    evecs.clear();
    double_vector_from_MatrixXcd(V, evecs);
    double_vector_from_VectorXcd(lambda, evals);
}

void eigen_decomposition_of_adjacency(graph& g, vector< std::complex<double> >& evals, vector< vector< std::complex<double> > >& evecs) {
    // Returns the Eigen Decomposition of the Graph Adjacency Matrix
    // Each Column of MatrixXcd represents 1 Eigen vector
    MatrixXd gmat;
    make_eigen_format_adjacency_matrix(g, gmat);
    VectorXcd lambda;
    MatrixXcd V;
    eigen_decomposition(gmat, lambda, V);
    evals.clear();
    evecs.clear();
    complex_vector_from_MatrixXcd(V, evecs);
    complex_vector_from_VectorXcd(lambda, evals);
}

void resistance_distance(graph& g, vector< vector<double> >& resistance) {
    // Returns the effective resistance between all pair of nodes, where the edge weight
    // represents the resistance between 2 given nodes.
    vector<double> evals;
    vector< vector<double> > evecs;
    eigen_decomposition_of_laplacian(g, evals, evecs);
    resistance.clear();
    resistance.assign(g.get_num_nodes(), vector<double> (g.get_num_nodes(), 0));
    for (unsigned long i = 0; i < g.get_num_nodes(); i++) {
        for (unsigned long j = 0; j < g.get_num_nodes(); j++) {
            for (unsigned long k = 0; k < evals.size(); k++) {
                if (evals[k] != 0) {
                    double t1 = evecs[i][k] - evecs[j][k];
                    resistance[i][j] += (1 / evals[k]) * (t1 * t1);
                }
            }
        }
    }
}

void conductance_distance(graph& g, vector< vector<double> >& conductance) {
    // Returns the effective conductance between all pair of nodes, where the edge weight
    // represents the conductance between 2 given nodes.
    resistance_distance(g, conductance);
    for (unsigned long i = 0; i < g.get_num_nodes(); i++)
        for (unsigned long j = 0; j < g.get_num_nodes(); j++)
            conductance[i][j] = 1 / conductance[i][j];
}

void hop_transition_probability(graph& g, vector< vector<double> >& transition_probability) {
    // Returns the Hop transition probability of a given graph as defined in
    // Paper : Unveiling Community Structures in Weighted Networks by Nelson A. Alves
    conductance_distance(g, transition_probability);
    for (unsigned long i = 0; i < g.get_num_nodes(); i++) {
        double denominator = 0;
        for (unsigned long j = 0; j < g.get_num_nodes(); j++) {
            if (i != j)
                denominator += transition_probability[i][j];
        }
        for (unsigned long j = 0; j < g.get_num_nodes(); j++)
            if (i != j)
                transition_probability[i][j] /= denominator;
    }
}

void average_distance(graph& g, vector< vector<double> >& hop_distance) {
    // Returns the average Hop distance between pair of nodes of a given graph as defined in
    // Paper : Unveiling Community Structures in Weighted Networks by Nelson A. Alves
    vector< vector<double> > transition_probability;
    hop_transition_probability(g, transition_probability);
    hop_distance.clear();
    hop_distance.assign(g.get_num_nodes(), vector<double> (g.get_num_nodes(), 0));
    for (unsigned long i = 0; i < g.get_num_nodes(); i++) {
        for (unsigned long j = 0; j < g.get_num_nodes(); j++) {
            double sum = 0;
            for (unsigned long k = 0; k < g.get_num_nodes(); k++) {
                if ((k != i) && (k != j))
                    sum += abs(transition_probability[i][k] - transition_probability[j][k]);
            }
            hop_distance[i][j] = sum / (g.get_num_nodes() - 2);
        }
    }
}

double efficiency_conductance(graph& g) {
    /* Returns the efficiency of the graph using pairwise node conductance.
     * Conductance is defined as the inverse of resistance. In graph theoretic sense,
     * the resistance between a pair of node is calculated as the effective resistance
     * in electrical network theory between the pair, when edge weights represents the
     * resistance between two nodes. Conductance Efficiency = sum of all pairwise
     * conductance of the graph/ sum of all pairwise conductance of a complete graph
     * of equivalent size. */
    graph complete(0, 0);
    double complete_conduct = 0, g_conduct = 0;
    generate_clique_graph(complete, g.get_num_nodes());
    vector< vector<double> > conductance;
    conductance_distance(complete, conductance);
    for (unsigned long i = 0; i < g.get_num_nodes(); i++) {
        for (unsigned long j = 0; j < g.get_num_nodes(); j++) {
            if (i != j)
                complete_conduct += conductance[i][j];
        }
    }
    conductance.clear();
    //    conductance_distance(g,conductance);
    resistance_distance(g, conductance);
    for (unsigned long i = 0; i < g.get_num_nodes(); i++) {
        for (unsigned long j = 0; j < g.get_num_nodes(); j++) {
            if (i != j)
                g_conduct += conductance[i][j];
        }
    }
    return (g_conduct / complete_conduct);
}


#endif	/* EIGEN_EXT_H */

