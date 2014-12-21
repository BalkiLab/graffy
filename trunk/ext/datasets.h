/*
 * File:   datasets.h
 * Author: sudip
 *
 * Created on 6 March, 2013, 12:02 AM
 */

#ifndef DATASETS_H
#define	DATASETS_H

#define ARTIFICIAL_DATASETS "/opt/Extensions/datasets/artificial/"
#define REAL_DATASETS "/opt/Extensions/datasets/real/"

void populate_p2p_graphs(vector<string>& files, bool ptype) {
    // ptype 1: fresh load all graphs to the vector   0: append all graphs to existing vector

    if (ptype) files.clear();
    files.push_back("chord_1000.el");
    files.push_back("chord_10000.el");
    files.push_back("db_1000.el");
    files.push_back("db_10000.el");
    files.push_back("kademlia_1000.el");
    files.push_back("kademlia_10000.el");
}

void populate_random_graphs(vector<string>& files, bool ptype) {
    // ptype 1: fresh load all graphs to the vector   0: append all graphs to existing vector

    if (ptype) files.clear();
    files.push_back("er_1000_0.1.el");
    files.push_back("er_5000_0.01.el");
}

void populate_small_world_graphs(vector<string>& files, bool ptype) {
    // ptype 1: fresh load all graphs to the vector   0: append all graphs to existing vector

    if (ptype) files.clear();
    files.push_back("sw_2_12_0.25.el");
    files.push_back("sw_3_8_0.25.el");
    files.push_back("sw_3_9_0.25.el");
}

void populate_scale_free_graphs(vector<string>& files, bool ptype) {
    // Generated from BA model of igraph
    // ptype 1: fresh load all graphs to the vector   0: append all graphs to existing vector

    if (ptype) files.clear();
    files.push_back("ba_2_1000.el");
    files.push_back("ba_2_10000.el");
    files.push_back("ba_3_1000.el");
    files.push_back("ba_3_10000.el");
}

void populate_forest_fire_graphs(vector<string>& files, bool ptype) {
    // ptype 1: fresh load all graphs to the vector   0: append all graphs to existing vector

    if (ptype) files.clear();
    files.push_back("ff_0.2_1000.el");
    files.push_back("ff_0.4_1000.el");
}

void populate_star_graphs(vector<string>& files, bool ptype) {
    // ptype 1: fresh load all graphs to the vector   0: append all graphs to existing vector

    if (ptype) files.clear();
    files.push_back("star_1000.el");
    files.push_back("star_5000.el");
    files.push_back("star_10000.el");
}

void populate_co_authors_networks(vector<string>& files, bool ptype) {
    // ptype 1: fresh load all graphs to the vector   0: append all graphs to existing vector

    if (ptype) files.clear();
    files.push_back("codmat99.el");
    files.push_back("codmat3.el");
    files.push_back("codmat5.el");
    files.push_back("astro-ph.el");
    files.push_back("ca-grqc.el");
    files.push_back("hep-th.el");
    files.push_back("netsci.el");
}

void populate_infrastructure_networks(vector<string>& files, bool ptype) {
    // ptype 1: fresh load all graphs to the vector   0: append all graphs to existing vector

    if (ptype) files.clear();
    files.push_back("power.el");
    files.push_back("asi.el");
    files.push_back("USairport_2010.el");
}

void populate_anti_social_networks(vector<string>& files, bool ptype) {
    // ptype 1: fresh load all graphs to the vector   0: append all graphs to existing vector

    if (ptype) files.clear();
    files.push_back("911.el");
    files.push_back("cni.el");
    files.push_back("gtt.el");
    files.push_back("2002.el");
    files.push_back("prison.el");
}

void populate_social_networks(vector<string>& files, bool ptype) {
    // ptype 1: fresh load all graphs to the vector   0: append all graphs to existing vector

    if (ptype) files.clear();
    files.push_back("emails.el");
    files.push_back("brightkite.el");
    files.push_back("gowalla.el");
}

void populate_team_networks(vector<string>& files, bool ptype) {
    // ptype 1: fresh load all graphs to the vector   0: append all graphs to existing vector

    if (ptype) files.clear();
    files.push_back("football.el");
    files.push_back("jazz.el");
}

void populate_purchase_networks(vector<string>& files, bool ptype) {
    // ptype 1: fresh load all graphs to the vector   0: append all graphs to existing vector

    if (ptype) files.clear();
    files.push_back("amazon.el");
    files.push_back("books.el");
}

void populate_ppi_networks(vector<string>& files, bool ptype) {
    // ptype 1: fresh load all graphs to the vector   0: append all graphs to existing vector

    if (ptype) files.clear();
    files.push_back("ppi_yeast.el");
    files.push_back("hprd.el");
}

#endif	/* DATASETS_H */

