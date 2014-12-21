/*
 * File:   R_ext.h
 * Author: bharath
 *
 * Created on 1 March, 2013, 10:17 PM
 */

#ifndef R_EXT_H
#define	R_EXT_H

#include "includes.h"

struct fitresult_R {
    double gof_tail_pl;
    id_type xmin;
    double alpha_tail_pl;
    double gof_all_pl;
    double alpha_all_pl;
    double gof_all_exp;
    double rate_all_exp;
    double gof_all_ln;
    double mean_all_ln;
    double sd_all_ln;
    double gof_head_pl;
    double alpha_head_pl;
    double gof_head_exp;
    double rate_head_exp;
    double gof_head_ln;
    double mean_head_ln;
    double sd_head_ln;

    fitresult_R() {
        gof_tail_pl = numeric_limits<double>::infinity();
        xmin = numeric_limits<long>::max();
        alpha_tail_pl = numeric_limits<double>::infinity();
        gof_all_pl = numeric_limits<double>::infinity();
        alpha_all_pl = numeric_limits<double>::infinity();
        gof_all_exp = numeric_limits<double>::infinity();
        rate_all_exp = numeric_limits<double>::infinity();
        gof_all_ln = numeric_limits<double>::infinity();
        mean_all_ln = numeric_limits<double>::infinity();
        sd_all_ln = numeric_limits<double>::infinity();
        gof_head_pl = numeric_limits<double>::infinity();
        alpha_head_pl = numeric_limits<double>::infinity();
        gof_head_exp = numeric_limits<double>::infinity();
        rate_head_exp = numeric_limits<double>::infinity();
        gof_head_ln = numeric_limits<double>::infinity();
        mean_head_ln = numeric_limits<double>::infinity();
        sd_head_ln = numeric_limits<double>::infinity();
    }

    void fit_degree_distribution(const graph& g, bool fit_indegrees) {
        char *temp_name;
        temp_name = tmpnam(NULL);
        string fileprefix(temp_name);
        ofstream ofs(fileprefix + ".in");
        for (id_type i = 0; i < g.get_num_nodes(); i++) {
            id_type deg_to_insert = (fit_indegrees) ? (double) g.get_node_in_degree(i) : (double) g.get_node_out_degree(i);
            if (deg_to_insert > 0)
                ofs << deg_to_insert << endl;
        }
        ofs.close();
        string command1 = "/opt/Extensions/R/bin/R CMD BATCH --no-save --no-restore '--args infilename=\"" + fileprefix + ".in\" outfilename=\"" + fileprefix + ".out\"' /opt/Extensions/RScripts/fitPowerlaw.R " + fileprefix + ".log";
        int code = system(command1.c_str());
        if (code == 0) {
            ifstream ifs(fileprefix + ".out");
            ifs >> gof_tail_pl;
            ifs >> xmin;
            ifs >> alpha_tail_pl;
            ifs >> gof_all_pl;
            ifs >> alpha_all_pl;
            ifs >> gof_all_exp;
            ifs >> rate_all_exp;
            ifs >> gof_all_ln;
            ifs >> mean_all_ln;
            ifs >> sd_all_ln;
            ifs >> gof_head_pl;
            ifs >> alpha_head_pl;
            ifs >> gof_head_exp;
            ifs >> rate_head_exp;
            ifs >> gof_head_ln;
            ifs >> mean_head_ln;
            ifs >> sd_head_ln;
            ifs.close();
        }
        string command2 = "rm " + fileprefix + ".in" + " " + fileprefix + ".out" + " " + fileprefix + ".log";
        code = system(command2.c_str());
    }
};

#endif	/* R_EXT_H */

