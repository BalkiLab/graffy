/* 
 * File:   graphio.h
 * Author: bharath
 *
 * Created on April 2, 2012, 11:46 AM
 */

#ifndef GRAPHIO_H
#define	GRAPHIO_H

#include <fstream>
#include <iostream>
#include <string>
#include "graph.h"
using namespace std;

namespace CDLib
{
        bool read_edgelist(graph& g,const string& filepath,bool directed, bool weighted);
        bool write_edgelist(graph& g,const string& filepath,bool weights);
        bool write_METIS(graph& g,const string& filepath,bool weights);
        bool write_SNAP(graph& g,const string& filepath,bool weights);
        bool write_SMAT(graph& g,const string& filepath,bool weights);
        bool write_UEL(graph& g,const string& filepath,bool weights);
        bool write_matlab_sp(graph& g,const string& filepath);
        bool write_dimacs_max_flow(graph& g, const string& filepath);
};

#endif	/* GRAPHIO_H */
