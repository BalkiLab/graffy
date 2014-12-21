/*
 * File:   graphio.cpp
 * Author: bharath
 *
 * Created on April 2, 2012, 11:46 AM
 */

#include "graphio.h"
#include "centrality.h"
using namespace CDLib;

bool CDLib::read_edgelist(graph& g, const string& filepath) {
    g.clear();
    ifstream ifs;
    ifs.open(filepath.c_str());
    if (ifs.is_open()) {
        vector<string> units;
        double weight;
        while (!ifs.eof()) {
            weight = 1;
            string line;
            getline(ifs, line);
            if ((line.size() > 0) && (line[0] != '#')) {
                split(line, units);
                if (units.size() != 0) {
                    if ((units.size() < 2) || (units.size() > 3)) return false;
                    if (units.size() == 3) weight = str2T<double>(units[2]);
                    g.add_node(units[0]);
                    g.add_node(units[1]);
                    g.add_edge(units[0], units[1], weight);
                }
            }
        }
        g.set_graph_name(filename(filepath));
        return true;
    }
    return false;
}

bool CDLib::read_adjacencylist(graph& g, const string& filepath) {
    g.clear();
    ifstream ifs;
    ifs.open(filepath.c_str());
    if (ifs.is_open()) {
        int type = 0;
        id_type nid = 0, estart = 0;
        string line;
        getline(ifs, line);
        vector<id_type> units;
        split(line, units);
        if ((units.size() > 2) && (units[0] == 0) && (units[1] == units.size() - 2)) {
            type = 0;
            estart = 2;
        } else if ((units.size() > 1) && (units[0] == 0)) {
            type = 1;
            estart = 1;
        } else {
            type = 2;
            estart = 0;
        }
        g.add_node(to_string(nid));
        for (id_type i = estart; i < units.size(); i++) {
            g.add_node(to_string(units[i]));
            g.add_edge(to_string(nid), to_string(units[i]), 1);
        }
        while (!ifs.eof()) {
            string line;
            getline(ifs, line);
            vector<id_type> units;
            split(line, units);
            if (units.size() > 0) {
                if ((type == 0) || (type == 1)) nid = units[0];
                else nid++;
                g.add_node(to_string(nid));
                for (id_type i = estart; i < units.size(); i++) {
                    g.add_node(to_string(units[i]));
                    g.add_edge(to_string(nid), to_string(units[i]), 1);
                }
            }
        }
        g.set_graph_name(filename(filepath));
        return true;
    }
    return false;
}

bool CDLib::read_matlab_sp(graph& g, const string& filepath) {
    g.clear();
    ifstream ifs;
    ifs.open(filepath.c_str());
    if (ifs.is_open()) {
        id_type from, to;
        double weight = 1;
        while (!ifs.eof()) {
            ifs >> from >> to >> weight;
            while (max(from, to) > g.get_num_nodes()) {
                g.add_node();
            }
            g.add_edge(from - 1, to - 1, weight);
        }
        g.set_graph_name(filename(filepath));
        return true;
    }
    return false;
}

bool CDLib::write_edgelist(const graph& g, const string& filepath, bool weights) {
    ofstream ofs;
    ofs.open(filepath.c_str());
    if (ofs.is_open()) {
        for (id_type i = 0; i < g.get_num_nodes(); i++) {
            for (adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i); aeit++) {
                ofs << g.get_node_label(i) << " " << g.get_node_label(aeit->first);
                if (weights) ofs << " " << aeit->second;
                ofs << endl;
            }
        }
        return true;
    }
    return false;
}

bool CDLib::write_xml(const graph& g, const string& filepath, bool weights) {
    ofstream ofs;
    string fp = filepath + ".xml";
    ofs.open(fp.c_str());
    if (ofs.is_open()) {
        string fid = filename(filepath);
        ofs << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>" << endl;
        ofs << "<DynamicMetaNetwork id=\"" << fid << "\">" << endl;
        ofs << "\t<MetaNetwork id=\"" << fid << "\">" << endl;
        ofs << "\t\t<nodes>" << endl;
        ofs << "\t\t\t<nodeclass type=\"Agent\" id=\"Agent\">" << endl;
        for (id_type i = 0; i < g.get_num_nodes(); i++) {
            ofs << "\t\t\t\t<node id=\"" << i << "\" title=\"" << g.get_node_label(i) << "\"/>" << endl;
        }
        ofs << "\t\t\t</nodeclass>" << endl;
        ofs << "\t\t</nodes>" << endl;
        ofs << "\t\t<networks>" << endl;
        ofs << "\t\t\t<graph sourceType=\"agent\" source=\"Agent\" targetType=\"agent\" target=\"Agent\" id=\"Agent x Agent\">" << endl;
        for (id_type i = 0; i < g.get_num_nodes(); i++) {
            for (adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i); aeit++) {
                ofs << "\t\t\t\t<edge source=\"" << i << "\" target=\"" << aeit->first << "\"";
                if (weights) ofs << " type=\"double\" value=\"" << aeit->second << "\"";
                ofs << "/>" << endl;
            }
        }
        ofs << "\t\t\t</graph>" << endl;
        ofs << "\t\t</networks>" << endl;
        ofs << "\t</MetaNetwork>" << endl;
        ofs << "</DynamicMetaNetwork>" << endl;
        return true;
    }
    return false;
}

bool CDLib::write_METIS(const graph& g, const string& filepath, bool weights) {
    ofstream file, file2;
    string graphfilename = filepath + ".metis";
    string labelfilename = filepath + ".metis_labels";
    file.open(graphfilename.c_str());
    file2.open(labelfilename.c_str());
    if (file.is_open() && file2.is_open()) {
        id_type edge_count = 0;
        for (id_type i = 0; i < g.get_num_nodes(); i++)
            for (adjacent_edges_iterator eit = g.out_edges_begin(i); eit != g.out_edges_end(i); eit++)
                if ((*eit).first != i) edge_count++;
        file << g.get_num_nodes() << " " << edge_count / 2 << " " << (weights ? "1" : "") << endl;
        for (id_type i = 0; i < g.get_num_nodes(); i++) {
            file2 << i << " " << g.get_node_label(i) << endl;
            for (adjacent_edges_iterator eit = g.out_edges_begin(i); eit != g.out_edges_end(i); eit++) {
                if ((*eit).first != i) {
                    file << ((*eit).first) + 1 << " ";
                    if (weights) file << (*eit).second << " ";
                }
            }
            file << endl;
        }
        file.close();
        return 1;
    } else return 0;
}

bool CDLib::write_SNAP(const graph& g, const string& filepath, bool weights) {
    ofstream file, file2;
    string graphfilename = filepath + ".snap";
    string labelfilename = filepath + ".snap_labels";
    file.open(graphfilename.c_str());
    file2.open(labelfilename.c_str());
    if (file.is_open() && file2.is_open()) {
        file << "p " << g.get_num_nodes() << " " << g.get_num_edges() << " " << ((g.is_directed()) ? "d" : "u") << " " << ((weights) ? "f" : "i") << " " << "0" << endl;
        for (id_type i = 0; i < g.get_num_nodes(); i++) {
            file2 << i << " " << g.get_node_label(i) << endl;
            for (adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i); aeit++)
                if (i <= aeit->first) file << i << " " << aeit->first << " " << aeit->second << endl;
        }
        return 1;
    } else return 0;
}

bool CDLib::write_SMAT(const graph& g, const string& filepath, bool weights) {
    ofstream file, file2;
    string graphfilename = filepath + ".smat";
    string labelfilename = filepath + ".smat_labels";
    file.open(graphfilename.c_str());
    file2.open(labelfilename.c_str());
    if (file.is_open() && file2.is_open()) {
        file << g.get_num_nodes() << " " << g.get_num_nodes() << " " << ((g.is_directed()) ? (g.get_num_edges() - g.get_num_self_edges()) : (2 * (g.get_num_edges() - g.get_num_self_edges()))) << endl;
        for (id_type i = 0; i < g.get_num_nodes(); i++) {
            file2 << i << " " << g.get_node_label(i) << endl;
            map<id_type, double> emap;
            for (adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i); aeit++)
                if (i != aeit->first)
                    emap.insert(*aeit);
            for (map<id_type, double>::iterator aeit = emap.begin(); aeit != emap.end(); aeit++)
                file << i << " " << aeit->first << " " << aeit->second << endl;
        }
        return 1;
    } else return 0;
}

bool CDLib::write_UEL(const graph& g, const string& filepath, bool weights) {
    ofstream file, file2;
    string graphfilename = filepath + ".net";
    string labelfilename = filepath + ".net_labels";
    file.open(graphfilename.c_str());
    file2.open(labelfilename.c_str());
    if (file.is_open() && file2.is_open()) {
        for (id_type i = 0; i < g.get_num_nodes(); i++) {
            file2 << i << " " << g.get_node_label(i) << endl;
            for (adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i); aeit++) {
                file << i << " " << aeit->first << " ";
                if (weights) file << " " << aeit->second;
                file << endl;
            }
        }
        return 1;
    } else return 0;
}

bool CDLib::write_matlab_sp(const graph& g, const string& filepath) {
    ofstream file;
    string graphfilename = filepath + ".spel";
    file.open(graphfilename.c_str());
    if (file.is_open()) {
        for (id_type i = 0; i < g.get_num_nodes(); i++) {
            for (adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i); aeit++)
                file << i + 1 << " " << aeit->first + 1 << " " << aeit->second << endl;
        }
        return 1;
    }
    return 0;
}

bool CDLib::write_dimacs_max_flow(const graph& g, const string& filepath) {
    ofstream file, file2;
    string graphfilename = filepath + ".dimacsflow";
    string labelfilename = filepath + ".dimacsflow_labels";
    file.open(graphfilename.c_str());
    file2.open(labelfilename.c_str());
    if (file.is_open() && file2.is_open()) {
        file << "p " << g.get_num_nodes() << " " << g.get_num_edges() << endl;
        file << "n 1 s" << endl;
        file << "n 2 t" << endl;
        for (id_type i = 0; i < g.get_num_nodes(); i++) {
            file2 << i + 1 << " " << g.get_node_label(i) << endl;
            for (adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i); aeit++)
                file << i + 1 << " " << aeit->first + 1 << " " << aeit->second << endl;
        }
        return 1;
    } else return 0;
}

bool CDLib::write_lcc_and_props(const graph& base, const string& filepath, bool start_with_one) {
    if (base.is_directed() || base.is_weighted()) return 0;
    node_set lcc;
    graph g(0, 0);
    cout << base.get_num_nodes() << " " << get_largest_connected_component(base, lcc) << endl;
    extract_subgraph(base, lcc, g);
    cout << base.get_num_edges() << " " << g.get_num_edges() << endl;
    ofstream file, file2, file3;
    string graphfilename = filepath + ".graph";
    string degdistfilename = filepath + ".degdist";
    string nccfilename = filepath + ".ncc";

    file.open(graphfilename.c_str());
    file2.open(degdistfilename.c_str());
    //file3.open(nccfilename.c_str());
    map<id_type, id_type> degdist;
    if (file.is_open() && file2.is_open())// && file3.is_open())
    {
        for (id_type i = 0; i < g.get_num_nodes(); i++) {
            set<id_type> neighbors;
            //file3 << node_clustering_coefficient(g,i) << endl;
            pair < map<id_type, id_type>::iterator, bool> ret = degdist.insert(make_pair(g.get_node_out_degree(i), 1));
            if (!ret.second)ret.first->second++;
            for (adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i); aeit++)
                neighbors.insert(aeit->first);
            for (set<id_type>::iterator it = neighbors.begin(); it != neighbors.end(); it++)
                file << i + start_with_one << " " << (*it) + start_with_one << endl;
        }
        for (map<id_type, id_type>::iterator it = degdist.begin(); it != degdist.end(); it++)
            file2 << it->first << " " << it->second << endl;
        return 1;
    } else return 0;
}
