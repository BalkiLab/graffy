/*
 * File:   octave_ext.h
 * Author: bharath
 *
 * Created on 1 March, 2013, 10:19 PM
 */

#ifndef OCTAVE_EXT_H
#define	OCTAVE_EXT_H

#include "includes.h"
namespace matlab_ext {
    const string MATLAB_INTERPRETER_PATH = "/usr/local/MATLAB/R2012b/bin/";
    int run_matlab_script(const string& script_filename) {
        string matlab_cmd = MATLAB_INTERPRETER_PATH + "matlab -nodisplay -nosplash -r try, " + script_filename + ", catch, exit(1), end, exit(0)";
        return system(matlab_cmd.c_str());
    }
}

string get_spectral(graph& g) {
    char *temp_name;
    temp_name = tmpnam(NULL);
    string tf(temp_name);
    string tmpfile = filename(tf);
    string script_dir = "/opt/Extensions/MScripts/";
    string matlab_exepath = "/usr/local/MATLAB/R2012b/bin/";
    string resultfile = tmpfile + "_result.txt";
    string filepath = tmpfile + "_graph";
    write_matlab_sp(g, filepath);
    string matlabfile = tmpfile + "_execute_matlab.m";
    string matlab_cmd = tmpfile + "_cmd.log";
    string command = "cp " + script_dir + "do_spectral.m " + "./";
    int status = system(command.c_str());
    if (status) perror("Problem coping do_spectral.m script to current directory.");
    ofstream mofs;
    mofs.open(matlabfile.c_str());
    if (mofs.is_open()) {
        mofs << "do_spectral('" + filepath + ".spel','" << g.get_graph_name() << "','" << resultfile << "');" << endl;
        mofs << "exit;" << endl;
    }
    mofs.close();
    command = matlab_exepath + "matlab -nodesktop -nosplash -nojvm -nodisplay < " + matlabfile + " >> " + matlab_cmd;
    status = system(command.c_str());
    if (status) perror("Problem executing matlab command line from C++ for this graph");
    string line;
    ifstream ifs;
    ifs.open(resultfile.c_str());
    getline(ifs,line);          // Remember the order of output "name spectral_radius spectral_gap largest_eval_laplacian"
    ifs.close();
    command = "rm do_spectral.m " + resultfile + " " + matlabfile + " " + filepath + ".spel " + filepath + ".spel_labels " + matlab_cmd;
    status = system(command.c_str());
    if (status) perror("Problem deleting temporary files");
    return line;
}

#endif	/* OCTAVE_EXT_H */

