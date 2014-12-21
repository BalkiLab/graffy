#ifndef LFR_H
#define	LFR_H

#include "includes.h"

    namespace lfr_impl {
        /*
         * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
         *                                                                               *
         *	This program is free software; you can redistribute it and/or modify     *
         *  it under the terms of the GNU General Public License as published by         *
         *  the Free Software Foundation; either version 2 of the License, or            *
         *  (at your option) any later version.                                          *
         *                                                                               *
         *  This program is distributed in the hope that it will be useful,              *
         *  but WITHOUT ANY WARRANTY; without even the implied warranty of               *
         *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
         *  GNU General Public License for more details.                                 *
         *                                                                               *
         *  You should have received a copy of the GNU General Public License            *
         *  along with this program; if not, write to the Free Software                  *
         *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA    *
         *                                                                               *
         * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
         *                                                                               *
         *  Created by Andrea Lancichinetti on 10/04/08 (email: arg.lanci@gmail.com)     *
         *	Modified on 30/03/09                                                     *
         *	Collaborators: Santo Fortunato and Filippo Radicchi                      *
         *  Location: ISI foundation, Turin, Italy                                       *
         *	Project: Benchmarking community detection programs                       *
         *                                                                               *
         * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
         */
        const int R2_IM1 = 2147483563;
        const int R2_IM2 = 2147483399;
        const int unlikely = -214741;
        const double R2_AM = (1.0 / R2_IM1);
        const int R2_IMM1 = (R2_IM1 - 1);
        const int R2_IA1 = 40014;
        const int R2_IA2 = 40692;
        const int R2_IQ1 = 53668;
        const int R2_IQ2 = 52774;
        const int R2_IR1 = 12211;
        const int R2_IR2 = 3791;
        const int R2_NTAB = 32;
        const int R2_NDIV = (1 + R2_IMM1 / R2_NTAB);
        const double R2_EPS = 1.2e-7;
        const double R2_RNMX = (1.0 - R2_EPS);
        long seed;

        void srand4(void) {
            seed = (long) time(NULL);
        }

        void srand5(int rank) {
            seed = (long) (rank);
        }

        double ran2(long *idum) {
            int j;
            long k;
            static long idum2 = 123456789;
            static long iy = 0;
            static long iv[R2_NTAB];
            double temp;
            if (*idum <= 0 || !iy) {
                if (-(*idum) < 1) *idum = 1 * (*idum);
                else *idum = -(*idum);
                idum2 = (*idum);
                for (j = R2_NTAB + 7; j >= 0; j--) {
                    k = (*idum) / R2_IQ1;
                    *idum = R2_IA1 * (*idum - k * R2_IQ1) - k*R2_IR1;
                    if (*idum < 0) *idum += R2_IM1;
                    if (j < R2_NTAB) iv[j] = *idum;
                }
                iy = iv[0];
            }
            k = (*idum) / R2_IQ1;
            *idum = R2_IA1 * (*idum - k * R2_IQ1) - k*R2_IR1;
            if (*idum < 0) *idum += R2_IM1;
            k = (idum2) / R2_IQ2;
            idum2 = R2_IA2 * (idum2 - k * R2_IQ2) - k*R2_IR2;
            if (idum2 < 0) idum2 += R2_IM2;
            j = iy / R2_NDIV;
            iy = iv[j] - idum2;
            iv[j] = *idum;
            if (iy < 1) iy += R2_IMM1;
            if ((temp = R2_AM * iy) > R2_RNMX) return R2_RNMX;
            else return temp;
        }

        double ran4(void) {
            return (ran2(&seed));
        }

        int irand(int n) {
            return (int(ran4()*(n + 1)));
        }

        void srand6(void) {
            srand4();
            seed = irand(R2_IM2);
            if (seed < 1 || seed > R2_IM2) seed = 1;
        }

        int powerlaw(int n, int min, double tau, deque<double> &cumulative) {
            cumulative.clear();
            double a = 0;
            for (double h = min; h < n + 1; h++)
                a += pow((1. / h), tau);
            double pf = 0;
            for (double i = min; i < n + 1; i++) {
                pf += 1 / a * pow((1. / (i)), tau);
                cumulative.push_back(pf);
            }
            return 0;
        }

        int config_model(deque <set<int> > &en, deque<int> &d, int min) {
            en.clear();
            set <int> first;
            for (int i = 0; i < d.size(); i++)
                en.push_back(first);
            multimap <int, int> degree_node;
            for (int i = 0; i < d.size(); i++)
                degree_node.insert(make_pair(d[i], i));
            int var = 0;
            while (degree_node.size() > 0) {
                multimap<int, int>::iterator itlast = degree_node.end();
                itlast--;
                multimap <int, int>::iterator itit = itlast;
                deque <multimap<int, int>::iterator> erasenda;
                int inserted = 0;
                for (int i = 0; i < itlast->first; i++) {
                    if (itit != degree_node.begin()) {
                        itit--;
                        en[itlast->second].insert(itit->second);
                        en[itit->second].insert(itlast->second);
                        inserted++;
                        erasenda.push_back(itit);
                    } else break;
                }
                for (int i = 0; i < erasenda.size(); i++) {
                    if (erasenda[i]->first > 1)
                        degree_node.insert(make_pair(erasenda[i]->first - 1, erasenda[i]->second));
                    degree_node.erase(erasenda[i]);
                }
                var += itlast->first - inserted;
                degree_node.erase(itlast);
            }
            for (int u = 0; u < en.size(); u++) {
                int stopper = 0;
                while (en[u].size() < min) {
                    if (stopper++>d.size() * d.size())return -1;
                    int nl = irand(en.size() - 1);
                    if (nl != u) {
                        en[u].insert(nl);
                        en[nl].insert(u);
                    }
                }
            }
            for (int nodo = 0; nodo < d.size(); nodo++) for (int krm = 0; krm < en[nodo].size(); krm++) {
                    int stopper2 = 0;
                    while (stopper2 < d.size()) {
                        stopper2++;
                        int old_node = 0;
                        int random_mate = irand(d.size() - 1);
                        while (random_mate == nodo)
                            random_mate = irand(d.size() - 1);
                        int nodo_h = 0;
                        if (en[nodo].insert(random_mate).second) {
                            deque <int> nodi_esterni;
                            for (set<int>::iterator it_est = en[nodo].begin(); it_est != en[nodo].end(); it_est++) if ((*it_est) != random_mate)
                                    nodi_esterni.push_back(*it_est);
                            old_node = nodi_esterni[irand(nodi_esterni.size() - 1)];
                            en[nodo].erase(old_node);
                            en[random_mate].insert(nodo);
                            en[old_node].erase(nodo);
                            deque <int> not_common;
                            for (set<int>::iterator it_est = en[random_mate].begin(); it_est != en[random_mate].end(); it_est++)
                                if ((old_node != (*it_est)) && (en[old_node].find(*it_est) == en[old_node].end()))
                                    not_common.push_back(*it_est);
                            nodo_h = not_common[irand(not_common.size() - 1)];
                            en[random_mate].erase(nodo_h);
                            en[nodo_h].erase(random_mate);
                            en[nodo_h].insert(old_node);
                            en[old_node].insert(nodo_h);
                            break;
                        }
                    }
                }
            return 0;
        }

        double integral(double a, double b) {
            if (fabs(a + 1.) > 1e-10) return (1. / (a + 1.) * pow(b, a + 1.));
            else return (log(b));
        }

        double average_degree(const double &dmax, const double &dmin, const double &gamma) {
            return (1. / (integral(gamma, dmax) - integral(gamma, dmin)))*(integral(gamma + 1, dmax) - integral(gamma + 1, dmin));
        }

        double solve_dmin(const double& dmax, const double &dmed, const double &gamma) {
            double dmin_l = 1;
            double dmin_r = dmax;
            double average_k1 = average_degree(dmin_r, dmin_l, gamma);
            double average_k2 = dmin_r;
            if ((average_k1 - dmed > 0) || (average_k2 - dmed < 0)) return -1;
            while (fabs(average_k1 - dmed) > 1e-7) {
                double temp = average_degree(dmax, ((dmin_r + dmin_l) / 2.), gamma);
                if ((temp - dmed)*(average_k2 - dmed) > 0) {
                    average_k2 = temp;
                    dmin_r = ((dmin_r + dmin_l) / 2.);
                } else {
                    average_k1 = temp;
                    dmin_l = ((dmin_r + dmin_l) / 2.);
                }
            }
            return dmin_l;
        }

        double integer_average(int n, int min, double tau) {
            double a = 0;
            for (double h = min; h < n + 1; h++)
                a += pow((1. / h), tau);
            double pf = 0;
            for (double i = min; i < n + 1; i++)
                pf += 1 / a * pow((1. / (i)), tau) * i;
            return pf;
        }

        int change_community_size(deque<int> &seq) {
            if (seq.size() <= 1)
                return 0;
            int min1 = 0;
            int min2 = 0;
            for (int i = 0; i < seq.size(); i++)
                if (seq[i] <= seq[min1])
                    min1 = i;
            if (min1 == 0)
                min2 = 1;
            for (int i = 0; i < seq.size(); i++)
                if (seq[i] <= seq[min2] && seq[i] > seq[min1])
                    min2 = i;
            seq[min1] += seq[min2];
            int c = seq[0];
            seq[0] = seq[min2];
            seq[min2] = c;
            seq.pop_front();
            return 0;
        }

        template<typename T, template<typename> class C>
        int shuffle_s(C<T> &sq) {
            int siz = sq.size();
            if (siz == 0)
                return -1;
            for (int i = 0; i < sq.size(); i++) {
                int random_pos = irand(siz - 1);
                T random_card_ = sq[random_pos];
                sq[random_pos] = sq[siz - 1];
                sq[siz - 1] = random_card_;
                siz--;
            }
            return 0;
        }

        int benchmark(CDLib::graph& g, vector<CDLib::id_type>& labels, int num_nodes, double average_k, int max_degree, double tau, double tau2, double mixing_parameter, int nmax = unlikely, int nmin = unlikely, bool excess = false, bool defect = false) {
            srand6();
            bool fixed_range = (nmax != unlikely);
            if (num_nodes == unlikely || average_k == unlikely || max_degree == unlikely || tau == unlikely || tau2 == unlikely || mixing_parameter == unlikely) return -1;
            if (num_nodes <= 0 || average_k <= 0 || max_degree <= 0 || mixing_parameter < 0 || (nmax <= 0 && nmax != unlikely) || (nmin <= 0 && nmin != unlikely)) return -1;
            if (fixed_range && (nmin > nmax || nmax > num_nodes)) return -1;
            if (excess && defect) return -1;
            if (mixing_parameter > 1) {
                excess = false;
                defect = false;
            }
            double dmin = solve_dmin(max_degree, average_k, -tau);
            if (dmin == -1) return -1;
            int min_degree = int(dmin);
            double media1 = integer_average(max_degree, min_degree, tau);
            double media2 = integer_average(max_degree, min_degree + 1, tau);
            if (fabs(media1 - average_k) > fabs(media2 - average_k)) min_degree++;
            if (!fixed_range) {
                nmax = max_degree;
                nmin = int(min_degree);
            }
            deque <int> degree_seq;
            deque <double> cumulative;
            powerlaw(max_degree, min_degree, tau, cumulative);
            for (int i = 0; i < num_nodes; i++) {
                int nn = lower_bound(cumulative.begin(), cumulative.end(), ran4()) - cumulative.begin() + min_degree;
                degree_seq.push_back(nn);
            }
            deque<set<int> > en;
            if (config_model(en, degree_seq, min_degree) == -1) return -1;
            for (int i = 0; i < en.size(); i++)
                degree_seq[i] = en[i].size();
            int max_degree_actual = 0;
            int max_degree_actual2 = 0;
            deque <int> internal_degree_seq;
            for (int i = 0; i < degree_seq.size(); i++) {
                double interno = (1 - mixing_parameter) * degree_seq[i];
                int int_interno = int(interno);
                if (ran4()<(interno - int_interno))
                    int_interno++;
                if (excess) {
                    while ((double(int_interno) / degree_seq[i] < (1 - mixing_parameter)) && (int_interno < degree_seq[i]))
                        int_interno++;
                }
                if (defect) {
                    while ((double(int_interno) / degree_seq[i] > (1 - mixing_parameter)) && (int_interno > 0))
                        int_interno--;
                }
                internal_degree_seq.push_back(int_interno);
                if (int_interno > max_degree_actual)
                    max_degree_actual = int_interno;
                if (degree_seq[i] > max_degree_actual2)
                    max_degree_actual2 = degree_seq[i];
            }
            deque <int> num_seq;
            powerlaw(nmax, nmin, tau2, cumulative);
            int _num_ = 0;
            if (!fixed_range && (max_degree_actual + 1) > nmin) {
                _num_ = max_degree_actual + 1;
                num_seq.push_back(max_degree_actual + 1);
            }
            while (true) {
                int nn = lower_bound(cumulative.begin(), cumulative.end(), ran4()) - cumulative.begin() + nmin;
                if (nn + _num_ <= num_nodes) {
                    num_seq.push_back(nn);
                    _num_ += nn;
                } else
                    break;
            }
            num_seq[num_seq.size() - 1] += num_nodes - _num_;
            int ncom = num_seq.size();
            deque<deque<int> > member_matrix;
            deque <int> first;
            for (int i = 0; i < ncom; i++) member_matrix.push_back(first);
            deque <int> refused(degree_seq.size());
            for (int i = 0; i < degree_seq.size(); i++) refused[i] = i;
            int k_r = 0;
            while (refused.size() > 0) {
                k_r++;
                deque<int> new_refused;
                for (int i = 0; i < refused.size(); i++) {
                    int random_module = irand(ncom - 1);
                    if (internal_degree_seq[refused[i]]<(num_seq[random_module])) {
                        if (member_matrix[random_module].size() == num_seq[random_module]) {
                            new_refused.push_back(member_matrix[random_module][0]);
                            member_matrix[random_module].pop_front();
                        }
                        member_matrix[random_module].push_back(refused[i]);
                    } else {
                        new_refused.push_back(refused[i]);
                    }
                }
                refused.clear();
                refused = new_refused;
                int missing_links = 0;
                for (int j = 0; j < refused.size(); j++)
                    missing_links += internal_degree_seq[refused[j]];
                if (k_r > 3 * num_nodes) {
                    k_r = 0;
                    ncom--;
                    member_matrix.pop_back();
                    for (int im = 0; im < member_matrix.size(); im++)
                        member_matrix[im].clear();
                    refused.clear();
                    for (int im = 0; im < degree_seq.size(); im++)
                        refused.push_back(im);
                    if (ncom == 0) return -1;
                    change_community_size(num_seq);
                }
            }

            if (ncom == 1 && mixing_parameter < 1) return -1;
            for (int i = 0; i < member_matrix.size(); i++) {
                int internal_cluster = 0;
                for (int j = 0; j < member_matrix[i].size(); j++)
                    internal_cluster += internal_degree_seq[member_matrix[i][j]];
                if (internal_cluster % 2 != 0) {
                    bool default_flag = false;
                    if (excess)
                        default_flag = true;
                    else if (defect)
                        default_flag = false;
                    else if (ran4() > 0.5)
                        default_flag = true;
                    if (default_flag) {
                        for (int j = 0; j < member_matrix[i].size(); j++) {
                            if ((internal_degree_seq[member_matrix[i][j]] < member_matrix[i].size() - 1) && (internal_degree_seq[member_matrix[i][j]] < degree_seq[member_matrix[i][j]])) {
                                internal_degree_seq[member_matrix[i][j]]++;
                                break;
                            }
                        }
                    } else {
                        for (int j = 0; j < member_matrix[i].size(); j++) {
                            if (internal_degree_seq[member_matrix[i][j]] > 0) {
                                internal_degree_seq[member_matrix[i][j]]--;
                                break;
                            }
                        }
                    }
                }
            }
            deque <int> member_list(num_nodes);
            for (int i = 0; i < member_matrix.size(); i++)
                for (int j = 0; j < member_matrix[i].size(); j++)
                    member_list[member_matrix[i][j]] = i;
            deque<int>internal_plusminus_seq(degree_seq.size());
            int var2 = 0;
            for (int i = 0; i < num_nodes; i++) {
                int inte = 0;
                for (set<int>::iterator it = en[i].begin(); it != en[i].end(); it++)
                    if (member_list[*it] == member_list[i])
                        inte++;
                internal_plusminus_seq[i] = inte - internal_degree_seq[i];
                var2 += (internal_plusminus_seq[i])*(internal_plusminus_seq[i]);
            }
            int best_var = var2;
            int max_counter = num_nodes / member_matrix.size();
            int counter_ = 0;
            if (mixing_parameter > 1) counter_ = max_counter;
            int deqar_size = 0;
            for (int i = 0; i < member_matrix.size(); i++)
                if (member_matrix[i].size() > deqar_size)
                    deqar_size = member_matrix[i].size();
            if (max_degree_actual2 > deqar_size)
                deqar_size = max_degree_actual2;
            deque<int> deqar(deqar_size);
            int deq_s;
            while (counter_ < max_counter) {
                counter_++;
                int stopper2_limit = 3;
                for (int i = 0; i < member_matrix.size(); i++) {
                    for (int j = 0; j < member_matrix[i].size(); j++) {
                        int nodo = member_matrix[i][j];
                        int stopper2 = 0;
                        while (stopper2 < stopper2_limit && internal_plusminus_seq[nodo] < 0) {
                            stopper2++;
                            deq_s = 0;
                            for (int k = 0; k < member_matrix[member_list[nodo]].size(); k++) {
                                int candidate = member_matrix[member_list[nodo]][k];
                                if ((en[nodo].find(candidate) == en[nodo].end()) && (candidate != nodo))
                                    deqar[deq_s++] = candidate;
                            }
                            if (deq_s == 0) break;
                            int random_mate = deqar[irand(deq_s - 1)];
                            {
                                deq_s = 0;
                                for (set<int>::iterator it_est = en[nodo].begin(); it_est != en[nodo].end(); it_est++)
                                    if (member_list[*it_est] != member_list[nodo])
                                        deqar[deq_s++] = *it_est;
                                if (deq_s == 0)
                                    break;
                                int old_node = deqar[irand(deq_s - 1)];
                                deq_s = 0;
                                for (set<int>::iterator it_est = en[random_mate].begin(); it_est != en[random_mate].end(); it_est++)
                                    if ((old_node != (*it_est)) && (en[old_node].find(*it_est) == en[old_node].end()))
                                        deqar[deq_s++] = *it_est;
                                if (deq_s == 0) break;
                                int nodo_h = deqar[irand(deq_s - 1)];
                                int nodes_here[4];
                                nodes_here[0] = nodo;
                                nodes_here[1] = random_mate;
                                nodes_here[2] = old_node;
                                nodes_here[3] = nodo_h;
                                int current_plus_minus[4];
                                for (int k = 0; k < 4; k++)
                                    current_plus_minus[k] = internal_plusminus_seq[nodes_here[k]];
                                current_plus_minus[0]++;
                                current_plus_minus[1]++;
                                if (member_list[nodo_h] == member_list[random_mate]) {
                                    current_plus_minus[1]--;
                                    current_plus_minus[3]--;
                                }
                                if (member_list[nodo_h] == member_list[old_node]) {
                                    current_plus_minus[2]++;
                                    current_plus_minus[3]++;
                                }
                                int impr = 0;
                                for (int k = 0; k < 4; k++)
                                    impr += (internal_plusminus_seq[nodes_here[k]])*(internal_plusminus_seq[nodes_here[k]]) - (current_plus_minus[k])*(current_plus_minus[k]);
                                if (impr < 0) break;
                                var2 -= impr;
                                en[nodo].erase(old_node);
                                en[nodo].insert(random_mate);
                                en[old_node].erase(nodo);
                                en[old_node].insert(nodo_h);
                                en[random_mate].erase(nodo_h);
                                en[random_mate].insert(nodo);
                                en[nodo_h].erase(random_mate);
                                en[nodo_h].insert(old_node);
                                for (int k = 0; k < 4; k++)
                                    internal_plusminus_seq[nodes_here[k]] = current_plus_minus[k];
                            }
                        }
                        while (stopper2 < stopper2_limit && internal_plusminus_seq[nodo] > 0) {
                            stopper2++;
                            int random_module = irand(member_matrix.size() - 1);
                            while (random_module == i)
                                random_module = irand(member_matrix.size() - 1);
                            int random_mate = member_matrix[random_module][irand(member_matrix[random_module].size() - 1)];
                            if (en[nodo].find(random_mate) == en[nodo].end()) {
                                deq_s = 0;
                                for (set<int>::iterator it_est = en[nodo].begin(); it_est != en[nodo].end(); it_est++)
                                    if (member_list[*it_est] == member_list[nodo])
                                        deqar[deq_s++] = *it_est;
                                if (deq_s == 0)
                                    break;
                                int old_node = deqar[irand(deq_s - 1)];
                                deq_s = 0;
                                for (set<int>::iterator it_est = en[random_mate].begin(); it_est != en[random_mate].end(); it_est++)
                                    if ((old_node != (*it_est)) && (en[old_node].find(*it_est) == en[old_node].end()))
                                        deqar[deq_s++] = *it_est;
                                if (deq_s == 0)
                                    break;
                                int nodo_h = deqar[irand(deq_s - 1)];
                                int nodes_here[4];
                                nodes_here[0] = nodo;
                                nodes_here[1] = random_mate;
                                nodes_here[2] = old_node;
                                nodes_here[3] = nodo_h;
                                int current_plus_minus[4];
                                for (int k = 0; k < 4; k++)
                                    current_plus_minus[k] = internal_plusminus_seq[nodes_here[k]];
                                current_plus_minus[0]--;
                                current_plus_minus[2]--;
                                if (member_list[nodo_h] == member_list[random_mate]) {
                                    current_plus_minus[1]--;
                                    current_plus_minus[3]--;
                                }
                                if (member_list[nodo_h] == member_list[old_node]) {
                                    current_plus_minus[2]++;
                                    current_plus_minus[3]++;
                                }
                                int impr = 0;
                                for (int k = 0; k < 4; k++)
                                    impr += (internal_plusminus_seq[nodes_here[k]])*(internal_plusminus_seq[nodes_here[k]]) - (current_plus_minus[k])*(current_plus_minus[k]);
                                if (impr < 0)
                                    break;
                                var2 -= impr;
                                en[nodo].erase(old_node);
                                en[nodo].insert(random_mate);
                                en[old_node].erase(nodo);
                                en[old_node].insert(nodo_h);
                                en[random_mate].erase(nodo_h);
                                en[random_mate].insert(nodo);
                                en[nodo_h].erase(random_mate);
                                en[nodo_h].insert(old_node);
                                for (int k = 0; k < 4; k++)
                                    internal_plusminus_seq[nodes_here[k]] = current_plus_minus[k];
                            }
                        }
                    }
                }
                if (var2 == 0)
                    break;
                if (best_var > var2) {
                    best_var = var2;
                    counter_ = 0;
                }
            }
            int eras_add_times = 0;
            if (excess) {
                for (int i = 0; i < num_nodes; i++)
                    while ((degree_seq[i] > 1) && double(internal_plusminus_seq[i] + internal_degree_seq[i]) / degree_seq[i] < 1 - mixing_parameter) {
                        int nodo = i;
                        ++eras_add_times;
                        deq_s = 0;
                        for (set<int>::iterator it_est = en[nodo].begin(); it_est != en[nodo].end(); it_est++)
                            if (member_list[*it_est] != member_list[nodo])
                                deqar[deq_s++] = *it_est;
                        if (deq_s == en[nodo].size())return -1;
                        int random_mate = deqar[irand(deq_s - 1)];
                        en[nodo].erase(random_mate);
                        en[random_mate].erase(nodo);
                        degree_seq[i]--;
                        degree_seq[random_mate]--;
                    }
            }
            if (defect) {
                for (int i = 0; i < num_nodes; i++)
                    while ((degree_seq[i] < num_nodes) && double(internal_plusminus_seq[i] + internal_degree_seq[i]) / degree_seq[i] > 1 - mixing_parameter) {
                        int nodo = i;
                        ++eras_add_times;
                        int stopper_here = num_nodes;
                        int stopper_ = 0;
                        int random_mate = irand(num_nodes - 1);
                        while (!(member_list[random_mate] != member_list[nodo] && en[nodo].find(random_mate) == en[nodo].end()) && (stopper_ < stopper_here)) {
                            random_mate = irand(num_nodes - 1);
                            stopper_++;
                        }
                        if (stopper_ == stopper_here) return -1;
                        en[nodo].insert(random_mate);
                        en[random_mate].insert(nodo);
                        degree_seq[i]++;
                        degree_seq[random_mate]++;
                    }
            }
            double ratio_true = 0;
            double variance = 0;
            deque <double> double_mixing;
            for (int i = 0; i < degree_seq.size(); i++) {
                double_mixing.push_back(1. - double(internal_plusminus_seq[i] + internal_degree_seq[i]) / degree_seq[i]);
                ratio_true += double(internal_plusminus_seq[i] + internal_degree_seq[i]) / degree_seq[i];
                variance += (double(internal_plusminus_seq[i] + internal_degree_seq[i]) / degree_seq[i])*(double(internal_plusminus_seq[i] + internal_degree_seq[i]) / degree_seq[i]);
            }
            ratio_true = ratio_true / degree_seq.size();
            variance = variance / degree_seq.size();
            variance -= ratio_true*ratio_true;
            double density = 0;
            double sparsity = 0;
            int edges = 0;
            for (int i = 0; i < member_matrix.size(); i++) {
                double media_int = 0;
                double media_est = 0;
                for (int j = 0; j < member_matrix[i].size(); j++) {
                    media_int += internal_plusminus_seq[member_matrix[i][j]] + internal_degree_seq[member_matrix[i][j]];
                    media_est += degree_seq[member_matrix[i][j]] - internal_plusminus_seq[member_matrix[i][j]] - internal_degree_seq[member_matrix[i][j]];
                    edges += degree_seq[member_matrix[i][j]];
                }
                double pair_num = (member_matrix[i].size()*(member_matrix[i].size() - 1));
                double pair_num_e = ((num_nodes - member_matrix[i].size())*(member_matrix[i].size()));
                density += media_int / pair_num;
                sparsity += media_est / pair_num_e;
            }
            density = density / member_matrix.size();
            sparsity = sparsity / member_matrix.size();
            g.clear();
            for (CDLib::id_type i = 0; i < num_nodes; i++) g.add_node();
            for (int u = 0; u < en.size(); u++) {
                set<int>::iterator it = en[u].begin();
                while (it != en[u].end())
                    g.add_edge(u, *(it++), 1);
            }
            g.remove_self_edges();
            labels.assign(g.get_num_nodes(), 0);
            for (int i = 0; i < member_list.size(); i++) labels[i] = member_list[i];
            return 0;
        }
    };
    namespace lfr {
        struct lfr_props{
            int num_nodes;
            double avg_degree;
            int max_degree;
            double degdist_exponent;
            double commsizedist_exponent;
            double mixing_parameter;
        };
        void lfr_model(const lfr_props& lfp, CDLib::graph& g, vector<CDLib::id_type>& labels) {
            g.clear();
            lfr_impl::benchmark(g, labels, lfp.num_nodes, lfp.avg_degree, lfp.max_degree, lfp.degdist_exponent, lfp.commsizedist_exponent, lfp.mixing_parameter, lfr::unlikely, lfr::unlikely, lfr::unlikely, lfr::unlikely);
        }
    }
#endif	/* LFR_H */

