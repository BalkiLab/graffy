/*
 * File:   svm_ext.h
 * Author: sudip
 *
 * Created on 3 March, 2013, 1:10 AM
 */

#ifndef SVM_EXT_H
#define	SVM_EXT_H

#include "includes.h"

struct svm_parameters {
    string svm_type, kernel_type;
    int polynomial_degree, number_of_classes, number_of_support_vectors;
    double gamma, coef;
    vector<double> rho, probA, probB;
    vector<int> label, nr_sv;
};

struct svm {
    string train_options, predict_options;
    string svm_path;
    int flag_pp, flag_model;
    struct svm_parameters model_parameters;
    double sigma;
    double mean_squared_error, squared_correlation_coefficient, accuracy;
    unsigned long correctly_classified, total_classified, cv_fold_number;
    double cv_mean_squared_error, cv_squared_correlation_coefficient, cv_accuracy;
    vector<double> predicted;
    vector< vector<double> > predicted_with_probability; // Valid when probability estimates for each class is requested in classification.

    svm() {
        svm_path = "/opt/Extensions/SVM/"; // Default SVM Path, svm-train and svm-predict are updated in this path to suite this function.
        train_options = "";
        predict_options = "";
        flag_pp = 0;
        flag_model = 0; // 0: C_SVC or NU_SVC     1: NU_SVR or EPSILON_SVR      2: one-class SVM
        sigma = 0; // Valid only when flag_pp == 1 and flag_model == 1
        mean_squared_error = 0; // Valid only when flag_model == 1
        squared_correlation_coefficient = 0; // Valid only when flag_model == 1
        accuracy = 0; // Valid only when flag_model == 0 or flag_model == 2
        correctly_classified = 0; // Valid only when flag_model == 0 or flag_model == 2
        total_classified = 0; // Valid only when flag_model == 0 or flag_model == 2
        cv_fold_number = 5; // This is used for cross validation accuracy (Default 5-fold cross validation)
        cv_mean_squared_error = 0; // Valid only when flag_model == 1
        cv_squared_correlation_coefficient = 0; // Valid only when flag_model == 1
        cv_accuracy = 0; // Valid only when flag_model == 0 or flag_model == 2
    }

    svm(map<char, double>& flags) {
        //        Pass a map for setting the flags of SVM train. Look into libsvm README for knowing the flags
        //        for eq, flag setting should be like flags["s"] = 0;      flags["t"] = 2;
        svm();
        ostringstream train, predict;
        for (map<char, double>::iterator it = flags.begin(); it != flags.end(); ++it) {
            if (it->first != 'q') {
                train << "-" << it->first << " " << it->second << " "; // Last space in mandatory.
                if (it->first == 's') {
                    if ((it->second == 3) || (it->second == 4))
                        flag_model = 1;
                    if (it->second == 2)
                        flag_model = 2;
                }
                if (it->first == 'b') {
                    predict << "b " << it->second << " "; // Last space in mandatory.
                    flag_pp = 1;
                }
                if (it->first == 'v') {
                    cv_fold_number = (unsigned int) it->second; // Setting cross validation number.
                }
            } else {
                train << "-" << it->first << " "; // Last space in mandatory.
            }
        }
        train_options = train.str();
        predict_options = predict.str();
    }

    void execute_svm(vector< vector<double> >& train_data, vector< vector<double> >& test_data) {
        char *temp_name;
        temp_name = tmpnam(NULL);
        string tmpfile(temp_name);
        int code;
        string train_file = tmpfile + ".traind";
        prepare_svm_train(train_data, train_file);
        string model_file = tmpfile + ".model";
        string test_file = tmpfile + ".testd";
        string stdout_stream = tmpfile + ".sout";
        string stderr_stream = tmpfile + ".serr";
        prepare_svm_train(test_data, test_file);
        string output_file = tmpfile + ".outf";
        string command0 = svm_path + "svm-train -v " + cv_fold_number + " -q " + train_file + " 1>" + stdout_stream + " 2>>" + stderr_stream;
        code = system(command0.c_str()); // Executing command0 i.e. Cross Validation
        if (code == 0) {
            get_cross_validation(stdout_stream);
        }
        // Always provide a model file name.
        string command1 = svm_path + "svm-train " + train_options + train_file + " " + model_file + " 1>" + stdout_stream + " 2>>" + stderr_stream;
        code = system(command1.c_str()); // Executing command1 i.e. Training
        if (code == 0) {
            get_model_parameters(model_file);
        }
        string command2 = svm_path + "svm-predict " + predict_options + test_file + " " + model_file + " " + output_file + " 1>" + stdout_stream + " 2>>" + stderr_stream;
        code = system(command2.c_str()); // Executing command2 i.e. Prediction
        if (code == 0) {
            get_predicted_labels(output_file);
            read_stdout_stream(stdout_stream);
        }
        string command3 = "rm " + train_file + " " + model_file + " " + test_file + " " + output_file + " " + stdout_stream + " " + stderr_stream;
        code = system(command3.c_str()); // Removing temporary files
    }

    void svm_train(vector< vector<double> >& train_data, string & model_file) {
        char *temp_name;
        temp_name = tmpnam(NULL);
        string tmpfile(temp_name);
        int code;
        string train_file = tmpfile + ".traind";
        prepare_svm_train(train_data, train_file);
        string stdout_stream = tmpfile + ".sout";
        string stderr_stream = tmpfile + ".serr";
        string command0 = svm_path + "svm-train -v " + cv_fold_number + " -q " + train_file + " 1>" + stdout_stream + " 2>>" + stderr_stream;
        code = system(command0.c_str()); // Executing command0 i.e. Cross Validation
        if (code == 0) {
            get_cross_validation(stdout_stream);
        }
        // Always provide a model file name.
        string command1 = svm_path + "svm-train " + train_options + train_file + " " + model_file + " 1>" + stdout_stream + " 2>>" + stderr_stream;
        code = system(command1.c_str()); // Executing command1 i.e. Training
        if (code == 0) {
            get_model_parameters(model_file);
        }
        string command3 = "rm " + train_file + " " + stdout_stream + " " + stderr_stream;
        code = system(command3.c_str()); // Removing temporary files
    }

    void svm_predict(vector< vector<double> >& test_data, string & model_file) {
        char *temp_name;
        temp_name = tmpnam(NULL);
        string tmpfile(temp_name);
        int code;
        string test_file = tmpfile + ".testd";
        string stdout_stream = tmpfile + ".sout";
        string stderr_stream = tmpfile + ".serr";
        prepare_svm_train(test_data, test_file);
        string output_file = tmpfile + ".outf";
        string command2 = svm_path + "svm-predict " + predict_options + test_file + " " + model_file + " " + output_file + " 1>" + stdout_stream + " 2>>" + stderr_stream;
        code = system(command2.c_str()); // Executing command2 i.e. Prediction
        if (code == 0) {
            get_predicted_labels(output_file);
            read_stdout_stream(stdout_stream);
        }
        string command3 = "rm " + test_file + " " + output_file + " " + stdout_stream + " " + stderr_stream;
        code = system(command3.c_str()); // Removing temporary files
    }
private:

    void prepare_svm_train(vector< vector<double> >& train_data, string train_file) {
        //        Each column represent a feature. The 0th column is the labels. For binary classification the labels should be +1 or -1. Rows are the data points.
        if (train_data.size() == 0)
            return;
        ofstream ofs(train_file.c_str()); // Train data temporary file in libsvm format.
        for (id_type i = 0; i < train_data.size(); i++) {
            ofs << (int) train_data[i][0]; // int casted for strictly using classification. real values there will lead to regression.
            for (id_type j = 1; j < train_data.size(); j++) {
                ofs << " " << j << ":" << train_data[i][j];
            }
            ofs << endl;
        }
        ofs.close();
    }

    void get_predicted_labels(string output_file) {
        ifstream ifs(output_file.c_str());
        if (ifs.is_open()) {
            double label;
            if ((flag_pp == 1) && (flag_model == 0)) {
                id_type num_classes;
                ifs >> num_classes;
                double val;
                while (!ifs.eof()) {
                    vector<double> data_point;
                    ifs >> label;
                    data_point.push_back(label);
                    for (id_type i = 0; i < num_classes; i++) {
                        ifs >> val;
                        data_point.push_back(val);
                    }
                    predicted.push_back(label);
                    predicted_with_probability.push_back(data_point);
                }
            } else {
                while (!ifs.eof()) {
                    ifs >> label;
                    predicted.push_back(label);
                }
            }
        }
        ifs.close();
    }

    void read_stdout_stream(string stdout_stream) {
        ifstream ifs(stdout_stream.c_str());
        if (ifs.is_open()) {
            if ((flag_pp == 1) && (flag_model == 1))
                ifs >> sigma;
            if (flag_model == 1)
                ifs >> mean_squared_error >> squared_correlation_coefficient;
            else
                ifs >> accuracy >> correctly_classified >> total_classified;
        }
        ifs.close();
    }

    void get_cross_validation(string stdout_stream) {
        ifstream ifs(stdout_stream.c_str());
        if (ifs.is_open()) {
            if (flag_model == 1)
                ifs >> cv_mean_squared_error >> cv_squared_correlation_coefficient;
            else
                ifs >> cv_accuracy;
        }
        ifs.close();
    }

    void get_model_parameters(string model_file) {
        ifstream ifs(model_file.c_str());
        if (ifs.is_open()) {
            string tmpdata;
            ifs >> tmpdata >> model_parameters.svm_type;
            ifs >> tmpdata >> model_parameters.kernel_type;
            if (model_parameters.kernel_type.compare("polynomial") == 0) {
                ifs >> tmpdata >> model_parameters.polynomial_degree;
            }
            if ((model_parameters.kernel_type.compare("polynomial") == 0) || (model_parameters.kernel_type.compare("rbf") == 0) || (model_parameters.kernel_type.compare("sigmoid") == 0)) {
                ifs >> tmpdata >> model_parameters.gamma;
            }
            if ((model_parameters.kernel_type.compare("polynomial") == 0) || (model_parameters.kernel_type.compare("sigmoid") == 0)) {
                ifs >> tmpdata >> model_parameters.coef;
            }
            ifs >> tmpdata >> model_parameters.number_of_classes;
            ifs >> tmpdata >> model_parameters.number_of_support_vectors;
            get_and_break_line<int>(ifs, model_parameters.label);
            get_and_break_line<double>(ifs, model_parameters.probA);
            get_and_break_line<double>(ifs, model_parameters.probB);
            get_and_break_line<int>(ifs, model_parameters.nr_sv);
        }
        ifs.close();
    }

    template <typename T>
    void get_and_break_line(istream& ifs, vector<T>& elems) {
        string s;
        getline(ifs, s);
        elems.clear();
        stringstream ss(s);
        string item;
        ss >> item;
        T value;
        while (!ss.eof()) {
            ss >> value;
            elems.push_back(value);
        }
    }

    string get_tmpfile_names(string & tmpfile) {
        string train_file = tmpfile + ".traind";
        string model_file = tmpfile + ".model";
        string test_file = tmpfile + ".testd";
        string stdout_stream = tmpfile + ".sout";
        string stderr_stream = tmpfile + ".serr";
        string output_file = tmpfile + ".outf";
        string created_files = "Train File: " + train_file + "\nModel File: " + model_file + "\nTest File: " + test_file;
        created_files += "\nSTDOUT Stream File: " + stdout_stream + "\nSTDERR Stream File: " + stderr_stream + "\nOutput File: " + output_file + "\n";
        return created_files;
    }
};

#endif	/* SVM_EXT_H */

