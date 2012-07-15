/* 
 * File:   graphio.cpp
 * Author: bharath
 * 
 * Created on April 2, 2012, 11:46 AM
 */

#include "graphio.h"
using namespace CDLib;

bool CDLib::read_edgelist(graph& g,const string& filepath,bool directed, bool weighted)
{
    ifstream ifs;
    ifs.open(filepath.c_str());
    if(ifs.is_open())
    {
        string label1,label2;
        double weight = 1;
        while(!ifs.eof())
        {
            ifs >> label1 >> label2;
            if(weighted) ifs >> weight;
            if(label1.size() && label2.size() && weight)
            {
                g.add_node(label1);
                g.add_node(label2);
                g.add_edge(label1,label2,weight);
                if(!directed)g.add_edge(label2,label1,weight);
            }
        }
        return true;
    }
    return false;
}

bool CDLib::write_edgelist(graph& g,const string& filepath,bool weights)
{
    ofstream ofs;
    ofs.open(filepath.c_str());
    if(ofs.is_open())
    {
        for(id_type i=0; i< g.get_num_nodes();i++)
        {
            for(adjacent_edges_iterator aeit = g.out_edges_begin(i);aeit != g.out_edges_end(i);aeit++)
            {
                ofs << g.get_node_label(i) << "\t" << g.get_node_label(aeit->first) ;
                if(weights) ofs << "\t" << aeit->second ;
                ofs << endl;
            }
        }
        return true;
    }
    return false;
}

bool CDLib::write_METIS(graph& g,const string& filepath,bool weights)
{
    ofstream file,file2;
    string graphfilename = filepath+ ".metis";
    string labelfilename = filepath+ ".metis_labels";
    file.open(graphfilename.c_str());
    file2.open(labelfilename.c_str());
    if(file.is_open() && file2.is_open())
    {
        id_type edge_count = 0;
        for(id_type i =0 ;i < g.get_num_nodes();i++)
            for(adjacent_edges_iterator eit = g.out_edges_begin(i); eit != g.out_edges_end(i);eit++)
                if((*eit).first != i) edge_count++;
                    file << g.get_num_nodes() << " " << edge_count/2<< " " << (weights ? "1":"") << endl; 
        for(id_type i =0 ;i < g.get_num_nodes();i++)
        {
            file2 << i <<" "<< g.get_node_label(i) << endl;
            for(adjacent_edges_iterator eit = g.out_edges_begin(i); eit != g.out_edges_end(i);eit++)
            {
                    if((*eit).first != i)
                    {
                        file << ((*eit).first) + 1<< " ";
                        if(weights) file << (*eit).second << " ";
                    }
            }
            file << endl;
        }
        file.close();
        return 1;
    }
    else return 0;
}

bool CDLib::write_SNAP(graph& g,const string& filepath,bool weights)
{
    ofstream file,file2;
    string graphfilename = filepath+ ".snap";
    string labelfilename = filepath+ ".snap_labels";
    file.open(graphfilename.c_str());
    file2.open(labelfilename.c_str());
    if(file.is_open() && file2.is_open())
    {
        file << "p " << g.get_num_nodes() << " " << g.get_num_edges() << " " << ((g.is_directed()) ? "d" : "u") << " " << ((weights) ? "f" : "i") << " " << "0" << endl;
        for(id_type i=0; i< g.get_num_nodes();i++)
        {
            file2 << i <<" "<< g.get_node_label(i) << endl;
            for(adjacent_edges_iterator aeit = g.out_edges_begin(i);aeit != g.out_edges_end(i);aeit++)
                if(i <= aeit->first) file <<i << " " << aeit->first << " " << aeit->second << endl;
        }
        return 1;
    }
    else return 0;
}

bool CDLib::write_matlab_sp(graph& g,const string& filepath)
{
    ofstream file,file2;
    string graphfilename = filepath+ ".spel";
    string labelfilename = filepath+ ".spel_labels";
    file.open(graphfilename.c_str());
    file2.open(labelfilename.c_str());
    if(file.is_open() && file2.is_open())
    {
        for(id_type i=0; i< g.get_num_nodes();i++)
        {
            file2 << i+1 <<" "<< g.get_node_label(i) << endl;
            for(adjacent_edges_iterator aeit = g.out_edges_begin(i);aeit != g.out_edges_end(i);aeit++)
                file <<i+1 << " " << aeit->first+1 << " " << aeit->second << endl;
        }
        return 1;
    }
    else return 0;    
}

bool CDLib::write_dimacs_max_flow(graph& g, const string& filepath)
{
    ofstream file,file2;
    string graphfilename = filepath+ ".dimacsflow";
    string labelfilename = filepath+ ".dimacsflow_labels";
    file.open(graphfilename.c_str());
    file2.open(labelfilename.c_str());
    if(file.is_open() && file2.is_open())
    {
        file << "p " << g.get_num_nodes() << " " << g.get_num_edges() << endl;
        file << "n 1 s" << endl;
        file << "n 2 t" << endl;
        for(id_type i=0; i< g.get_num_nodes();i++)
        {
            file2 << i+1 <<" "<< g.get_node_label(i) << endl;
            for(adjacent_edges_iterator aeit = g.out_edges_begin(i);aeit != g.out_edges_end(i);aeit++)
                file <<i+1 << " " << aeit->first+1 << " " << aeit->second << endl;
        }
        return 1;
    }
    else return 0; 
}
