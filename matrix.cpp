/* 
 * File:   matrix.cpp
 * Author: sudip
 *
 * Created on 22 May, 2012, 10:21 PM
 */

#define EIGEN_USE_MKL_ALL
#include "matrix.h"
#include "graph.h"
#include <eigen3/Eigen/Dense>

using namespace CDLib;
using namespace Eigen;

void make_eigen_format_adjacency_matrix(graph& g, MatrixXd& gmat)
{
    // Creates the adjacency matrix for the given graph
    gmat.resize(g.get_num_nodes(),g.get_num_nodes());
    if (g.is_directed())
    {
        for(id_type i=0; i<g.get_num_nodes(); i++)
            for(id_type j=0; j<g.get_num_nodes(); j++)
                gmat(i,j) =  g.get_edge_weight(i,j);
    }
    else
    {
        for(id_type i=0; i<g.get_num_nodes(); i++){
            for(id_type j=0; j<g.get_num_nodes(); j++){
                gmat(i,j) =  g.get_edge_weight(i,j);
                gmat(j,i) =  g.get_edge_weight(j,i);
            }
        }
    }
}
 
void make_eigen_format_laplacian_matrix(graph& g, MatrixXd& lmat)
{
    // Creates the laplacian matrix for the given graph
    lmat.resize(g.get_num_nodes(),g.get_num_nodes());
    if (g.is_directed())
    {
        for(id_type i=0; i<g.get_num_nodes(); i++){
            for(id_type j=0; j<g.get_num_nodes(); j++){
                if(g.get_edge_weight(i,j))
                    lmat(i,j) = -1;
                else
                    lmat(i,j) = 0;
            }
        }
    }
    else
    {
        for(id_type i=0; i<g.get_num_nodes(); i++){
            for(id_type j=0; j<g.get_num_nodes(); j++){
                if(g.get_edge_weight(i,j)){
                    lmat(i,j) = -1;
                    lmat(j,i) = -1;
                }
                else{
                    lmat(i,j) = 0;
                    lmat(j,i) = 0;
                }
            }
        }
    }
    for(id_type i=0; i<g.get_num_nodes(); i++)
        lmat(i,i) = g.get_node_out_degree(i);
}

void eigen_decomposition(MatrixXd& mat, VectorXcd& lambda, MatrixXcd& V)
{
    // Creates the eigen decomposition of the given matrix into eigen values and eigen vectors 
    if (mat.rows() == mat.cols())
    {
        EigenSolver<MatrixXd> es(mat);
        lambda = es.eigenvalues();
        V = es.eigenvectors();
    }          
}

void complex_vector_from_VectorXcd(VectorXcd& v,vector< std::complex<double> >& vecs)
{
    // Convert VectorXcd type in Eigen Library to Vector of Complex Double type of C++ STL
    long size = v.size();
    for (long i=0; i<size; i++)
        vecs.push_back(v(i));
}

void complex_vector_from_MatrixXcd(MatrixXcd& v,vector< vector< std::complex<double> > >& mats)
{
    // Convert MatrixXcd type in Eigen Library to Vector of Vector of Complex  Double type of C++ STL
    // Each Column represents 1 Eigen vector
    long rows = v.rows();
    for (long i=0; i<rows; i++)
    {
        vector< std::complex<double> > temp;
        VectorXcd operation_vs = v.row(i);
        complex_vector_from_VectorXcd(operation_vs,temp);
        mats.push_back(temp);
    }
}

void double_vector_from_VectorXcd(VectorXcd& v,vector<double>& vecs)
{
    // Convert VectorXcd type in Eigen Library to Vector of Double type of C++ STL
    long size = v.size();
    for (long i=0; i<size; i++)
        vecs.push_back(real(v(i)));
}
 
void double_vector_from_MatrixXcd(MatrixXcd& v,vector< vector<double> >& mats)
{
    // Convert MatrixXcd type in Eigen Library to Vector of Vector of Double type of C++ STL
    // Each Column represents 1 Eigen vector
    long rows = v.rows();
    for (long i=0; i<rows; i++)
    {
        vector<double> temp;
        VectorXcd operation_vs = v.row(i);
        double_vector_from_VectorXcd(operation_vs,temp);
        mats.push_back(temp);
    }
}
 
void CDLib::eigen_decomposition_of_laplacian(graph& g,vector<double>& evals, vector< vector<double> >& evecs)
{
    // Returns the Eigen Decomposition of the Graph Laplacian
    MatrixXd lmat;
    make_eigen_format_laplacian_matrix(g,lmat);
    VectorXcd lambda;
    MatrixXcd V;
    eigen_decomposition(lmat,lambda,V);
    evals.clear();
    evecs.clear();
    double_vector_from_MatrixXcd(V,evecs);
    double_vector_from_VectorXcd(lambda,evals);
}
 
void CDLib::eigen_decomposition_of_adjacency(graph& g, vector< std::complex<double> >& evals, vector< vector< std::complex<double> > >& evecs)
{
    // Returns the Eigen Decomposition of the Graph Adjacency Matrix
    // Each Column of MatrixXcd represents 1 Eigen vector
    MatrixXd gmat;
    make_eigen_format_adjacency_matrix(g,gmat);
    VectorXcd lambda;
    MatrixXcd V;
    eigen_decomposition(gmat,lambda,V);
    evals.clear();
    evecs.clear();
    complex_vector_from_MatrixXcd(V,evecs);
    complex_vector_from_VectorXcd(lambda,evals);
}

void CDLib::resistance_distance(graph& g, vector< vector<double> >& resistance)
{
    // Returns the effective resistance between all pair of nodes, where the edge weight 
    // represents the resistance between 2 given nodes.
    vector<double> evals;
    vector< vector<double> > evecs;
    eigen_decomposition_of_laplacian(g,evals,evecs);
    resistance.clear();
    resistance.assign(g.get_num_nodes(),vector<double> (g.get_num_nodes(),0));
    for (unsigned long i=0; i<g.get_num_nodes(); i++){
        for (unsigned long j=0; j<g.get_num_nodes(); j++){
            for (unsigned long k=0; k<evals.size(); k++){
                if (evals[k] != 0){
                    double t1 = evecs[i][k] - evecs[j][k];
                    resistance[i][j] += (1/evals[k]) * (t1*t1);
                }
            }
        }
    }
}

void CDLib::conductance_distance(graph& g, vector< vector<double> >& conductance)
{
    // Returns the effective conductance between all pair of nodes, where the edge weight 
    // represents the conductance between 2 given nodes.
    resistance_distance(g,conductance);
    for (unsigned long i=0; i<g.get_num_nodes(); i++)
        for (unsigned long j=0; j<g.get_num_nodes(); j++)
            conductance[i][j] = 1/conductance[i][j];
}

void CDLib::hop_transition_probability(graph& g, vector< vector<double> >& transition_probability)
{
    // Returns the Hop transition probability of a given graph as defined in 
    // Paper : Unveiling Community Structures in Weighted Networks by Nelson A. Alves
    conductance_distance(g,transition_probability);
    for (unsigned long i=0; i<g.get_num_nodes(); i++){
        double denominator = 0;
        for (unsigned long j=0; j<g.get_num_nodes(); j++){
            if (i!=j)
                denominator += transition_probability[i][j];
        }
        for (unsigned long j=0; j<g.get_num_nodes(); j++)
            if (i!=j)
                transition_probability[i][j] /= denominator;
    }
}

void CDLib::average_distance(graph& g, vector< vector<double> >& hop_distance)
{
    // Returns the average Hop distance between pair of nodes of a given graph as defined in 
    // Paper : Unveiling Community Structures in Weighted Networks by Nelson A. Alves
    vector< vector<double> > transition_probability;
    hop_transition_probability(g,transition_probability);
    hop_distance.clear();
    hop_distance.assign(g.get_num_nodes(),vector<double> (g.get_num_nodes(),0));
    for (unsigned long i=0; i<g.get_num_nodes(); i++){
        for (unsigned long j=0; j<g.get_num_nodes(); j++){
            double sum = 0;
            for (unsigned long k=0; k<g.get_num_nodes(); k++){
                if ((k!=i) && (k!=j))
                    sum += abs(transition_probability[i][k] - transition_probability[j][k]);
            }
            hop_distance[i][j] = sum/(g.get_num_nodes()-2);
        }
    }
}

