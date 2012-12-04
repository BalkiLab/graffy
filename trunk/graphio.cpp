/* 
 * File:   graphio.cpp
 * Author: bharath
 * 
 * Created on April 2, 2012, 11:46 AM
 */

#include "graphio.h"
#include "centrality.h"
using namespace CDLib;

bool CDLib::read_edgelist(graph& g,const string& filepath,bool directed, bool weighted)
{
    if(directed!=g.is_directed() && weighted!=g.is_weighted()) return false;
    g.clear();
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

bool CDLib::read_matlab_sp(graph& g,const string& filepath){
    g.clear();
    ifstream ifs;
    ifs.open(filepath.c_str());
    if(ifs.is_open())
    {
        id_type from,to;
        double weight = 1;
        while(!ifs.eof())
        {
            ifs >> from >> to >> weight; 
            while( max(from,to) > g.get_num_nodes()){
                g.add_node();
            }
            g.add_edge(from-1,to-1,weight);
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
                ofs << g.get_node_label(i) << " " << g.get_node_label(aeit->first) ;
                if(weights) ofs << " " << aeit->second ;
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

bool CDLib::write_SMAT(graph& g,const string& filepath,bool weights)
{
    ofstream file,file2;
    string graphfilename = filepath+ ".smat";
    string labelfilename = filepath+ ".smat_labels";
    file.open(graphfilename.c_str());
    file2.open(labelfilename.c_str());
    if(file.is_open() && file2.is_open())
    {
        file << g.get_num_nodes() << " " << g.get_num_nodes() << " " << ((g.is_directed())?(g.get_num_edges()-g.get_num_self_edges()) : (2*(g.get_num_edges()-g.get_num_self_edges()))) << endl;
        for(id_type i=0; i< g.get_num_nodes();i++)
        {
            file2 << i <<" "<< g.get_node_label(i) << endl;
            map<id_type,double> emap;
            for(adjacent_edges_iterator aeit = g.out_edges_begin(i);aeit != g.out_edges_end(i);aeit++)
                if(i != aeit->first)
                    emap.insert(*aeit);
            for(map<id_type,double>::iterator aeit = emap.begin(); aeit != emap.end();aeit++)
                 file <<i << " " << aeit->first<< " " << aeit->second << endl;
        }
        return 1;
    }
    else return 0;    
}

bool CDLib::write_UEL(graph& g,const string& filepath,bool weights)
{
    ofstream file,file2;
    string graphfilename = filepath+ ".net";
    string labelfilename = filepath+ ".net_labels";
    file.open(graphfilename.c_str());
    file2.open(labelfilename.c_str());
    if(file.is_open() && file2.is_open())
    {
        for(id_type i=0; i< g.get_num_nodes();i++)
        {
            file2 << i <<" "<< g.get_node_label(i) << endl;
            for(adjacent_edges_iterator aeit = g.out_edges_begin(i);aeit != g.out_edges_end(i);aeit++)
            {
                file <<i << " " << aeit->first << " " ;
                if(weights) file << " " << aeit->second ;
                file << endl;
            }
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

bool CDLib::write_lcc_and_props(graph& base,const string& filepath,bool start_with_one){
    if(base.is_directed() || base.is_weighted()) return 0;
    node_set lcc;
    graph g(0,0);
    cout <<base.get_num_nodes()<< " "<< get_largest_connected_component(base,lcc) << endl;
    extract_subgraph(base,lcc,g);
    cout <<base.get_num_edges()<< " "<< g.get_num_edges() << endl;
    ofstream file,file2,file3;
    string graphfilename = filepath+ ".graph";
    string degdistfilename = filepath+ ".degdist";
    string nccfilename = filepath+ ".ncc"; 
    
    file.open(graphfilename.c_str());
    file2.open(degdistfilename.c_str());
    //file3.open(nccfilename.c_str());
    map<id_type,id_type> degdist;
    if(file.is_open() && file2.is_open())// && file3.is_open())
    {
        for(id_type i=0; i< g.get_num_nodes();i++)
        {
            set<id_type> neighbors;
            //file3 << node_clustering_coefficient(g,i) << endl;
            pair<map<id_type,id_type>::iterator,bool> ret = degdist.insert(make_pair(g.get_node_out_degree(i),1));
            if(!ret.second)ret.first->second++;
            for(adjacent_edges_iterator aeit = g.out_edges_begin(i);aeit != g.out_edges_end(i);aeit++)
                neighbors.insert(aeit->first);
            for(set<id_type>::iterator it = neighbors.begin();it!=neighbors.end();it++)
                file << i+start_with_one << " " << (*it) + start_with_one <<endl;
        }
        for(map<id_type,id_type>::iterator it = degdist.begin(); it != degdist.end();it++)
            file2 << it->first << " " << it->second << endl;
        return 1;
    }
    else return 0; 
}