/*
 * File:   utitlity.h
 * Author: bharath
 *
 * Created on May 13, 2012, 6:04 AM
 */

#ifndef UTITLITY_H
#define	UTITLITY_H

#include "typedefs.h"
#include "random.h"
using namespace CDLib;

template <typename T>
T str2T(const string& str) {
    T retval;
    istringstream iss(str);
    iss >> retval;
    return retval;
}

template <typename T>
string T2str(T val) {
    ostringstream oss;
    oss << val;
    return oss.str();
}

template <typename T>
void split(const string& s, const char delim, vector<T>& elems) {
    elems.clear();
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        if (item.length()) elems.push_back(str2T<T > (item));
        item.clear();
    }
}

template <typename T>
void split(const string& s, vector<T>& elems) {
    elems.clear();
    stringstream ss(s);
    string item;
    while (!ss.eof()) {
        ss >> item;
        if (item.length()) elems.push_back(str2T<T > (item));
        item.clear();
    }
}

template <typename T>
void cumsum(const vector<T>& v_in, vector<T>& v_out) {
    v_out.assign(v_in.size(), 0);
    v_out[0] = v_in[0];
    for (CDLib::id_type i = 1; i < v_in.size(); i++)
        v_out[i] = v_out[i - 1] + v_in[i];
}

template <class Container>
id_type num_unique_elements_seqential(const Container &c) {
    unordered_set<typename Container::value_type> test;
    for (typename Container::const_iterator it = c.begin(); it != c.end(); it++)
        test.insert(*it);
    return test.size();
}

template<typename T>
void read_vector_1D(const string& infile, vector<T>& data) {
    ifstream ifs(infile);
    copy(istream_iterator<T > (ifs), istream_iterator<T > (), back_inserter(data));
    ifs.close();
}

template<typename T>
void read_vector_of_vector(const string& infile, vector< vector<T> >& data) {
    ifstream ifs(infile);
    while (!ifs.eof()) {
        string line;
        getline(ifs, line);
        if(line != "" && line[0] != ' ' && line[0] != '\t' && line[0] != '\n'){
            istringstream iss(line);
            data.push_back(vector<T > ());
            while (!iss.eof()) {
                T val;
                iss >> val;
                data.back().push_back(val);
            }
        }
    }
    ifs.close();
}

template <class Element, class Val>
class max_label_picker {
private:
    unordered_map<Element, Val> label_fitness_map;
    vector<Element> max_labels;
    CDLib::id_type max_label_index;
    Val max_val;
    unordered_set<Element> curr_lbls;
public:
    typedef typename vector<Element>::const_iterator max_labels_iterator;
    typedef typename unordered_map<Element, Val>::const_iterator all_labels_iterator;

    max_label_picker(id_type node_degree) {
        max_val = (numeric_limits<Val>::has_infinity) ? -numeric_limits<Val>::infinity() : -numeric_limits<Val>::max();
        max_label_index = 0;
        max_labels.assign(node_degree, ((numeric_limits<Element>::has_infinity) ? -numeric_limits<Element>::infinity() : -numeric_limits<Element>::max()));
    }

    inline id_type size() const {
        return max_label_index;
    }

    inline bool empty() const {
        return max_label_index == 0;
    }

    inline id_type num_labels() const {
        return label_fitness_map.size();
    }

    inline Val max_value() const {
        return max_val;
    }

    void insert(pair<id_type, double> lbl_fit) {
        pair < unordered_map<id_type, double>::iterator, bool> ret = label_fitness_map.insert(lbl_fit);
        if (!ret.second)ret.first->second += lbl_fit.second;
        if (ret.first->second >= max_val) {
            if (ret.first->second > max_val) {
                max_label_index = 0;
                curr_lbls.clear();
                max_val = ret.first->second;
            }
            if (!max_label_index || curr_lbls.find(ret.first->first) == curr_lbls.end()) {
                max_labels[max_label_index] = ret.first->first;
                max_label_index++;
                curr_lbls.insert(ret.first->first);
            }
        }
    }

    inline typename vector<Element>::const_iterator max_labels_begin() const {
        if (max_label_index)return max_labels.begin();
        else return max_labels.end();
    }

    inline typename vector<Element>::const_iterator max_labels_end() const {
        if (max_label_index)return max_labels.begin() + max_label_index;
        else return max_labels.end();
    }

    inline typename vector<Element>::const_iterator find_in_max_labels(Element elem) const {
        if (max_label_index) return find(max_labels.begin(), max_labels.begin() + max_label_index, elem);
        else return max_labels.end();
    }

    inline typename vector<Element>::const_iterator get_random_max_label() {
        if (max_label_index != 0) {
            RandomGenerator<id_type> gen(0, max_label_index, 0);
            return max_labels.begin() + gen.next();
        } else return max_labels.end();
    }

    inline typename unordered_map<Element, Val>::const_iterator all_labels_begin() {
        return label_fitness_map.begin();
    }

    inline typename unordered_map<Element, Val>::const_iterator all_labels_end() {
        return label_fitness_map.end();
    }

    inline typename unordered_map<Element, Val>::const_iterator find_in_all_labels(id_type label) {
        return label_fitness_map.find(label);
    }
};

template <class Element, typename Priority>
struct min_heap {

    bool operator()(const pair<Element, Priority>& left, const pair<Element, Priority>& right) const {
        return (left.second < right.second);
    }
};

template <class Element, typename Priority>
struct max_heap {

    bool operator()(const pair<Element, Priority>& left, const pair<Element, Priority>& right) const {
        return (left.second > right.second);
    }
};

template <class Element, class Priority, class Compare>
class hash_heap {
private:
    deque< pair<Element, Priority> >vpiw_heap;
    unordered_map<Element, id_type> umis_pos;
    Compare cmp;
public:
    typedef typename deque< pair<Element, Priority> >::const_iterator iterator;

    id_type size() const {
        return vpiw_heap.size();
    }

    bool empty() const {
        return !vpiw_heap.size();
    }

    pair<Element, Priority> front() const {
        return vpiw_heap[0];
    }

    iterator begin() const {
        return vpiw_heap.begin();
    }

    iterator end() const {
        return vpiw_heap.end();
    }

    iterator find(Element elem) const {
        typename unordered_map<Element, id_type>::const_iterator it = umis_pos.find(elem);
        if (it != umis_pos.end()) return vpiw_heap.begin() + it->second;
        else return end();
    }

    bool update_key(const pair<Element, Priority>& piw_in) {
        if (umis_pos.find(piw_in.first) != umis_pos.end()) {
            id_type pos = umis_pos[piw_in.first];
            pair<Element, Priority> old_val = vpiw_heap[pos];
            vpiw_heap[pos].second = piw_in.second;
            if (!cmp(piw_in, old_val))heapify_up(pos);
            else heapify_down(pos);
            return true;
        }
        return insert(piw_in);
    }

    bool insert(const pair<Element, wt_type>& piw_in) {
        if (umis_pos.find(piw_in.first) != umis_pos.end()) return update_key(piw_in);
        vpiw_heap.push_back(piw_in);
        id_type pos = vpiw_heap.size() - 1;
        umis_pos.insert(make_pair(piw_in.first, pos));
        update_key(piw_in);
        return true;
    }

    void heapify_up(id_type pos) {
        while (pos > 0 && cmp(vpiw_heap[pos], vpiw_heap[(pos - 1) / 2])) {
            swap(umis_pos[vpiw_heap[pos].first], umis_pos[vpiw_heap[(pos - 1) / 2].first]);
            swap(vpiw_heap[pos], vpiw_heap[(pos - 1) / 2]);
            pos = (pos - 1) / 2;
        }
    }

    void heapify_down(id_type root) {
        if (root < vpiw_heap.size()) {
            id_type left = 2 * root + 1, right = 2 * root + 2, top = root;
            if (left < vpiw_heap.size() && cmp(vpiw_heap[left], vpiw_heap[top])) top = left;
            if (right < vpiw_heap.size() && cmp(vpiw_heap[right], vpiw_heap[top])) top = right;
            if (top != root) {
                swap(umis_pos[vpiw_heap[top].first], umis_pos[vpiw_heap[root].first]);
                swap(vpiw_heap[top], vpiw_heap[root]);
                heapify_down(top);
            }
        }
    }

    void pop() {
        umis_pos.erase(vpiw_heap[0].first);
        umis_pos[vpiw_heap[vpiw_heap.size() - 1].first] = 0;
        vpiw_heap[0] = vpiw_heap[vpiw_heap.size() - 1];
        vpiw_heap.pop_back();
        heapify_down(0);
    }

    bool remove(Element elem) {
        if (umis_pos.find(elem) == umis_pos.end()) return false;
        id_type pos = umis_pos[elem];
        umis_pos.erase(vpiw_heap[pos].first);
        umis_pos[vpiw_heap[vpiw_heap.size() - 1].first] = pos;
        vpiw_heap[pos] = vpiw_heap[vpiw_heap.size() - 1];
        vpiw_heap.pop_back();
        heapify_down(pos);
    }

    hash_heap() : vpiw_heap(), umis_pos(), cmp() {
    }

    hash_heap(const vector< pair<Element, Priority> >& v) {
        for (id_type i = 0; i < v.size(); i++) {
            vpiw_heap.push_back(make_pair(v[i].first, v[i].second));
            umis_pos.insert(make_pair(v[i].first, i));
        }
        for (id_type i = (vpiw_heap.size() + 1) / 2; i > 0; i--) heapify_down(i);
        heapify_down(0);
    }
};

typedef hash_heap<id_type, double, min_heap<id_type, double> > id_dbl_min_heap;
typedef hash_heap<id_type, double, max_heap<id_type, double> > id_dbl_max_heap;

template <typename K, typename V>
inline void map_insert_and_increment(unordered_map<K, V>& data, const K& key, const V& value) {
    pair < typename unordered_map<K, V>::iterator, bool> ret = data.insert(make_pair(key, 0));
    ret.first->second += value;
}

template <typename K, typename V>
inline V map_find_value_or_zero(const unordered_map<K, V>& data, const K& key) {
    typename unordered_map<K, V>::const_iterator ret = data.find(key);
    if (ret != data.end()) return ret->second;
    return 0;
}

template <typename K, typename V>
inline void map_find_and_modify(unordered_map<K, V>& data, const K& key, const V& value) {
    typename unordered_map<K, V>::iterator it = data.find(key);
    if (it != data.end()) it->second += value;
}

template <typename K, typename V>
inline void map_find_and_modify_force(unordered_map<K, V>& data, const K& key, const V& value) {
    typename unordered_map<K, V>::iterator it = data.find(key);
    if (it != data.end()) it->second += value;
    else map_insert_and_increment<K, V > (data, key, value);
}

inline string filename(string file) {
    ostringstream oss;
    size_t found_dot, found_slash;
    found_dot = file.find_last_of(".");
    found_slash = file.find_last_of("/");
    if ((found_dot != string::npos) && (found_slash != string::npos))
        oss << file.substr(found_slash + 1, (found_dot - found_slash - 1));
    else if ((found_dot == string::npos) && (found_slash != string::npos))
        oss << file.substr(found_slash + 1, (file.size() - found_slash - 1));
    else if ((found_dot != string::npos) && (found_slash == string::npos))
        oss << file.substr(0, found_dot);
    else
        oss << file;
    return oss.str();
}

inline void progress(int i, string comments) {
    cout << "\rProcessed " << i << " " << comments;
    fflush(stdout);
}

inline string change_extension(string file, string ext) {
    ostringstream oss;
    size_t found;
    found = file.find_last_of(".");
    if (found != string::npos)
        oss << file.substr(0, found) << "." << ext;
    else
        oss << file << "." << ext;
    return oss.str();
}

inline string get_extension(string file) {
    ostringstream oss;
    size_t found;
    found = file.find_last_of(".");
    if (found != string::npos)
        oss << file.substr(found + 1, string::npos);
    else
        oss << "";
    return oss.str();
}

inline string print_1d_vector(const vector<double>& mats, const string separator) {
    ostringstream oss;
    oss << endl;
    for (id_type i = 0; i < mats.size(); i++) {
        oss << fixed << mats[i] << separator;
    }
    return oss.str();
}

inline string print_2d_vector(const vector< vector<double> >& mats) {
    ostringstream oss;
    oss << endl;
    for (id_type i = 0; i < mats.size(); i++) {
        for (id_type j = 0; j < mats[i].size(); j++)
            oss << fixed << mats[i][j] << "\t";
        oss << endl;
    }
    return oss.str();
}

template <typename T1, typename T2>
bool compare_second_asc(const pair<T1, T2>& i, const pair<T1, T2>& j) {
    return i.second < j.second;
}

template <typename T1, typename T2>
bool compare_second_desc(const pair<T1, T2>& i, const pair<T1, T2>& j) {
    return i.second > j.second;
}

inline void sorting(vector< pair<id_type, double> >& unified, bool asc) {
    //    Sorts the pair in ascending order of the second value when bool asc is set else descending.
    if (asc)
        sort(unified.begin(), unified.end(), compare_second_asc<id_type, double>);
    else
        sort(unified.begin(), unified.end(), compare_second_desc<id_type, double>);
}

inline void delete_files(string location, vector<string>& files) {
    for (unsigned int i = 0; i < files.size(); i++) {
        string rem = location + files[i];
        remove(rem.c_str());
    }
}

#endif	/* UTITLITY_H */

