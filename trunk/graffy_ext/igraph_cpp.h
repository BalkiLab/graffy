/* 
 * File:   igraph_cpp.h
 * Author: bharath
 *
 * Created on 4 March, 2013, 1:23 AM
 */

#ifndef IGRAPH_CPP_H
#define	IGRAPH_CPP_H

#include "includes.h"
#include "igraph.h"
#include "igraph_ext.h"

namespace igcpp {

    class igcpp_igraph_t {
    private:
        igraph_t _igraph_t_internal;
        string graph_name;

        int import_graph(const graph& g) {
            igraph_vector_t edges = IGRAPH_VECTOR_NULL;
            IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
            IGRAPH_CHECK(igraph_vector_reserve(&edges, g.get_num_edges()));
            for (id_type i = 0; i < g.get_num_nodes(); i++) {
                for (adjacent_edges_iterator aeit = g.out_edges_begin(i); aeit != g.out_edges_end(i); aeit++) {
                    if (i <= aeit->first) {
                        IGRAPH_CHECK(igraph_vector_push_back(&edges, static_cast<long> (i)));
                        IGRAPH_CHECK(igraph_vector_push_back(&edges, static_cast<long> (aeit->first)));
                    }
                }
            }
            IGRAPH_CHECK(igraph_create(&_igraph_t_internal, &edges, g.get_num_nodes(), g.is_directed()));
            graph_name = g.get_graph_name();
            igraph_vector_destroy(&edges);
            IGRAPH_FINALLY_CLEAN(1);
            return 1;
        }

        int export_graph(CDLib::graph& g) const {
            g.clear();
            for (int k = 0; k < igraph_vcount(&_igraph_t_internal); k++)g.add_node();
            igraph_eit_t it;
            IGRAPH_CHECK(igraph_eit_create(&_igraph_t_internal, igraph_ess_all(IGRAPH_EDGEORDER_FROM), &it));
            IGRAPH_FINALLY(igraph_eit_destroy, &it);
            while (!IGRAPH_EIT_END(it)) {
                igraph_integer_t from, to;
                igraph_edge(&_igraph_t_internal, IGRAPH_EIT_GET(it), &from, &to);
                g.add_edge(static_cast<id_type> (from), static_cast<id_type> (to), 1);
                IGRAPH_EIT_NEXT(it);
            }
            g.set_graph_name(graph_name);
            igraph_eit_destroy(&it);
            IGRAPH_FINALLY_CLEAN(1);
            return 1;
        }

    public:

        //Unimplemented Functions
        //int igraph_degree(const igraph_t *graph, igraph_vector_t *res, const igraph_vs_t vids igraph_neimode_t mode, igraph_bool_t loops);
        //int igraph_add_edges(igraph_t *graph, const igraph_vector_t *edges, void *attr);
        //int igraph_add_vertices(igraph_t *graph, igraph_integer_t nv, void *attr);
        //int igraph_delete_edges(igraph_t *graph, igraph_es_t edges);
        //int igraph_delete_vertices(igraph_t *graph, const igraph_vs_t vertices);
        //int igraph_incident(const igraph_t *graph, igraph_vector_t *eids, igraph_integer_t pnode, igraph_neimode_t mode);
        //int igraph_neighbors(const igraph_t *graph, igraph_vector_t *neis, igraph_integer_t pnode, igraph_neimode_t mode);
        //int igraph_get_eids_multi(const igraph_t *graph, igraph_vector_t *eids,const igraph_vector_t *pairs, const igraph_vector_t *path,igraph_bool_t directed, igraph_bool_t error);
        //int igraph_get_eids(const igraph_t *graph, igraph_vector_t *eids, const igraph_vector_t *pairs, const igraph_vector_t *path, igraph_bool_t directed, igraph_bool_t error);
        //int igraph_get_eid(const igraph_t *graph, igraph_integer_t *eid, igraph_integer_t pfrom, igraph_integer_t pto, igraph_bool_t directed, igraph_bool_t error);
        //int igraph_edge(const igraph_t *graph, igraph_integer_t eid, igraph_integer_t *from, igraph_integer_t *to);


        //To pass to igraph functions

        inline igraph_t* c_igraph_t() {
            return &_igraph_t_internal;
        }

        inline const igraph_t* c_igraph_t() const {
            return &_igraph_t_internal;
        }

        //igraph_integer_t igraph_vcount(const igraph_t *graph);

        inline id_type get_num_nodes() const {
            return igraph_vcount(&_igraph_t_internal);
        }

        //igraph_integer_t igraph_ecount(const igraph_t *graph);

        inline id_type get_num_edges() const {
            return igraph_ecount(&_igraph_t_internal);
        }

        //igraph_bool_t igraph_is_directed(const igraph_t *graph);

        inline bool is_directed() const {
            return igraph_is_directed(&_igraph_t_internal);
        }

        inline void set_graph_name(const string& name) {
            graph_name = name;
        }

        inline string& get_graph_name(const string& name) {
            return graph_name;
        }

        //To a graffy graph

        inline bool graffy_graph(graph& g) const {
            return export_graph(g);
        }


        //Construct from graffy_graph

        inline igcpp_igraph_t(const graph& g) {
            import_graph(g);
        }

        //int igraph_destroy(igraph_t *graph);

        inline ~igcpp_igraph_t() {
            igraph_destroy(&_igraph_t_internal);
        }

        //int igraph_empty(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed); Empty Undirected Graph

        inline igcpp_igraph_t() {
            igraph_empty(&_igraph_t_internal, 0, IGRAPH_UNDIRECTED);
        }

        //int igraph_empty(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed);

        inline igcpp_igraph_t(id_type num_nodes, bool directed) {
            igraph_empty(&_igraph_t_internal, num_nodes, ((directed) ? IGRAPH_DIRECTED : IGRAPH_UNDIRECTED));
        }

        //Copy Constructor int igraph_copy(igraph_t *to, const igraph_t *from);

        inline igcpp_igraph_t(const igcpp_igraph_t& ig2) {
            igraph_copy(&_igraph_t_internal, ig2.c_igraph_t());
        }
    };

    //const igraph_vector_t*igraph_vector_view (const igraph_vector_t *v, const igraph_real_t *data, long int length);

//    inline igraph_vector_t* to_igraph_vector_t(vector<igraph_real_t>& v) {
//        igraph_vector_t iv;
//        return igraph_vector_view(&iv, &v[0], v.size());
//    }

    class igcpp_igraph_vector_t {
    private:
        igraph_vector_t _internal_igraph_vector_t;
    public:
        typedef igraph_real_t* iterator;
        typedef const igraph_real_t* const_iterator;

        //void igraph_vector_clear(igraph_vector_t* v);

        void clear() {
            igraph_vector_clear(&_internal_igraph_vector_t);
        }

        //int igraph_vector_insert(igraph_vector_t *v, long int pos,igraph_real_t value);

        void insert(id_type pos, igraph_real_t val) {
            igraph_vector_insert(&_internal_igraph_vector_t, pos, val);
        }

        void insert(const_iterator it, igraph_real_t val) {
            igraph_vector_insert(&_internal_igraph_vector_t, it - begin(), val);
        }

        //void igraph_vector_remove(igraph_vector_t *v, long int elem);

        void erase(const_iterator it) {
            igraph_vector_remove(&_internal_igraph_vector_t, it - begin());
        }

        void erase(id_type pos) {
            igraph_vector_remove(&_internal_igraph_vector_t, pos);
        }

        void erase(id_type start_pos, id_type end_pos) {
            igraph_vector_remove_section(&_internal_igraph_vector_t, start_pos, end_pos);
        }

        void erase(const_iterator it1, const_iterator it2) {
            igraph_vector_remove_section(&_internal_igraph_vector_t, static_cast<id_type>(it1-begin()), static_cast<id_type>(it2-begin()));
        }

        //igraph_real_t igraph_vector_pop_back(igraph_vector_t* v);

        void pop_back() {
            igraph_vector_pop_back(&_internal_igraph_vector_t);
        }

        //int igraph_vector_push_back (igraph_vector_t* v, igraph_real_t e);

        void push_back(igraph_real_t val) {
            igraph_vector_push_back(&_internal_igraph_vector_t, val);
        }

        igraph_real_t& operator[] (id_type n) {
            return VECTOR(_internal_igraph_vector_t)[n];
        }

        const igraph_real_t& operator[] (id_type n) const {
            return VECTOR(_internal_igraph_vector_t)[n];
        }

        //igraph_real_t* igraph_vector_e_ptr  (const igraph_vector_t* v, long int pos);

        iterator begin() {
            return igraph_vector_e_ptr(&_internal_igraph_vector_t, 0);
        }

        iterator end() {
            return igraph_vector_e_ptr(&_internal_igraph_vector_t, size() - 1);
        }

        //igraph_real_t* igraph_vector_e_ptr  (const igraph_vector_t* v, long int pos);

        const_iterator begin() const {
            return igraph_vector_e_ptr(&_internal_igraph_vector_t, 0);
        }

        const_iterator end() const {
            return igraph_vector_e_ptr(&_internal_igraph_vector_t, size() - 1);
        }


        //To pass to igraph functions

        inline igraph_vector_t* c_igraph_vector_t() {
            return &_internal_igraph_vector_t;
        }

        inline const igraph_vector_t* c_igraph_vector_t() const {
            return &_internal_igraph_vector_t;
        }
        
        inline bool empty() const {
            return (size() == 0);
        }

        //long int igraph_vector_size      (const igraph_vector_t* v);

        inline id_type size() const {
            return igraph_vector_size(&_internal_igraph_vector_t);
        }

        //long int igraph_vector_capacity(const igraph_vector_t*v);

        inline id_type capacity() const {
            return igraph_vector_capacity(&_internal_igraph_vector_t);
        }

        //void igraph_vector_fill(igraph_vector_t* v, igraph_real_t e);

        inline void fill(igraph_real_t val) {
            igraph_vector_fill(&_internal_igraph_vector_t, val);
        }

        //int igraph_vector_resize(igraph_vector_t* v, long int newsize);

        inline void resize(id_type size) {
            igraph_vector_resize(&_internal_igraph_vector_t, size);
        }

        //int igraph_vector_reserve   (igraph_vector_t* v, long int size);

        inline void reserve(id_type size) {
            igraph_vector_reserve(&_internal_igraph_vector_t, size);
        }

        inline void assign(id_type size, igraph_real_t val) {
            igraph_vector_fill(&_internal_igraph_vector_t, val);
        }

        //void igraph_vector_fill (igraph_vector_t* v, igraph_real_t e);

        inline igcpp_igraph_vector_t(){
            igraph_vector_init(&_internal_igraph_vector_t, 0);
        }
        
        inline igcpp_igraph_vector_t(id_type size, igraph_real_t val) {
            igraph_vector_init(&_internal_igraph_vector_t, 0);
            assign(size, val);
        }

        //int igraph_vector_copy(igraph_vector_t *to, const igraph_vector_t *from); - Copy Constructor

        inline igcpp_igraph_vector_t(const igcpp_igraph_vector_t& v) {
            igraph_vector_init(&_internal_igraph_vector_t, 0);
            igraph_vector_copy(&_internal_igraph_vector_t, v.c_igraph_vector_t());
        }

        //int igraph_vector_init_copy(igraph_vector_t *v, igraph_real_t *data, long int length);

        inline igcpp_igraph_vector_t(const vector<igraph_real_t>& v) {
            igraph_vector_init(&_internal_igraph_vector_t, v.size());
            for(id_type i=0;i<v.size();i++) VECTOR(_internal_igraph_vector_t)[i] = v[i];
        }

        template<typename T>
        inline void stl_vector(vector<T>&v) const{
            v.assign(size(), 0);
            for (id_type i = 0; i < size(); i++)
                v[i] = static_cast<T> (VECTOR(_internal_igraph_vector_t)[i]);

        }

        //int igraph_vector_init_seq(igraph_vector_t *v, igraph_real_t from, igraph_real_t to);

        inline igcpp_igraph_vector_t(igraph_real_t start, igraph_real_t end) {
            igraph_vector_init_seq(&_internal_igraph_vector_t, start, end);
        }

        // void igraph_vector_destroy(igraph_vector_t* v); - Destructor

        inline ~igcpp_igraph_vector_t() {
            igraph_vector_destroy(&_internal_igraph_vector_t);
        }
    };

    class igcpp_igraph_matrix_t {
    private:
        igraph_matrix_t _internal_igraph_matrix_t;
    public:

        inline igraph_matrix_t* c_igraph_matrix_t(){
            return &_internal_igraph_matrix_t;
        }
        
        inline const igraph_matrix_t* c_igraph_matrix_t() const{
            return &_internal_igraph_matrix_t;
        }
        
        ~igcpp_igraph_matrix_t(){
            igraph_matrix_destroy(&_internal_igraph_matrix_t);
        }
        
        inline void assign(long int nrow, long int ncol, igraph_real_t e) {
            igraph_matrix_fill(&_internal_igraph_matrix_t,e);
        }

        inline igcpp_igraph_matrix_t() {
            igraph_matrix_init(&_internal_igraph_matrix_t, 0, 0);
        }

        inline igcpp_igraph_matrix_t(long int nrow, long int ncol, igraph_real_t e) {
            igraph_matrix_init(&_internal_igraph_matrix_t, 0, 0);
            assign(nrow, ncol, e);
        }

        inline igcpp_igraph_matrix_t(const igcpp_igraph_matrix_t& m) {
            igraph_matrix_init(&_internal_igraph_matrix_t, 0, 0);
            igraph_matrix_copy(&_internal_igraph_matrix_t, m.c_igraph_matrix_t());
        }

        inline id_type nrow() const{
            return igraph_matrix_nrow(&_internal_igraph_matrix_t);
        }
        
        inline id_type ncol() const{
            return igraph_matrix_ncol(&_internal_igraph_matrix_t);
        }
        
        inline double& get_value(id_type row, id_type col){
            return *igraph_matrix_e_ptr(&_internal_igraph_matrix_t,row, col);
        }

    };


}

#endif	/* IGRAPH_CPP_H */

