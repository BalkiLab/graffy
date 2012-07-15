#ifndef EXTERN_HPP
#define EXTERN_HPP
#include "typedefs.h"
namespace lfr_binary_networks 
{    
/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                               *
*	This program is free software; you can redistribute it and/or modify         *
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
*  Created by Andrea Lancichinetti on 7/01/09 (email: arg.lanci@gmail.com)      *
*	Modified on 28/05/09                                                         *
*	Collaborators: Santo Fortunato												 *
*  Location: ISI foundation, Turin, Italy                                       *
*	Project: Benchmarking community detection programs                           *
*                                                                               *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/
    #define R2_IM1 2147483563
    #define R2_IM2 2147483399
    #define R2_AM (1.0/R2_IM1)
    #define R2_IMM1 (R2_IM1-1)
    #define R2_IA1 40014
    #define R2_IA2 40692
    #define R2_IQ1 53668
    #define R2_IQ2 52774
    #define R2_IR1 12211
    #define R2_IR2 3791
    #define R2_NTAB 32
    #define R2_NDIV (1+R2_IMM1/R2_NTAB)
    #define R2_EPS 1.2e-7
    #define R2_RNMX (1.0-R2_EPS)
    #define unlikely -214741

    long cherr() {
        cerr << "the check failed" << endl;
        long e;
        cin >> e;
        return e;
    }
    long cherr(double a) {
        cerr << "the check failed because of " << a << endl;
        long e;
        cin >> e;
        return e;
    }
    template <typename uno, typename due>
    void prints(pair <uno, due> &sq, ostream &out) {
        out << sq.first << "\t" << sq.second << endl;
    }
    template <typename uno, typename due>
    void prints(pair <uno, due> &sq) {
        cout << sq.first << "\t" << sq.second << endl;
    }
    template <typename uno, typename due>
    void prints(map <uno, due> &sq, ostream &out) {
        typename map <uno, due>::iterator it = sq.begin();
        while (it != sq.end()) {
            out << it->first << "\t" << it->second << endl;
            it++;
        }
        out << endl;
    }
    template <typename uno, typename due>
    void prints(multimap <uno, due> &sq, ostream &out) {
        typename map <uno, due>::iterator it = sq.begin();
        while (it != sq.end()) {
            out << it->first << "\t" << it->second << endl;
            it++;
        }
        out << endl;
    }
    template <typename Seq>
    void prints(Seq &sq, ostream &out) {
        typename Seq::iterator it = sq.begin();
        while (it != sq.end())
            out << *(it++) << "\t";
        out << endl;
    }
    template <typename type_>
    void prints(type_ *a, long b) {
        for (long i = 0; i < b; i++)
            cout << a[i] << " ";
        cout << endl;
    }
    template<typename T, template<typename> class C>
    void printm(C<T>& c, ostream &out) {
        typename C<T>::iterator it = c.begin();
        while (it != c.end()) {
            prints(*it, out);
            it++;
        }
        out << endl;
    }
    template <typename uno, typename due>
    void prints(map <uno, due> &sq) {
        typename map <uno, due>::iterator it = sq.begin();
        while (it != sq.end()) {
            cout << it->first << "\t" << it->second << endl;
            it++;
        }
        cout << endl;
    }
    template <typename uno, typename due>
    void prints(multimap <uno, due> &sq) {
        typename map <uno, due>::iterator it = sq.begin();
        while (it != sq.end()) {
            cout << it->first << "\t" << it->second << endl;
            it++;
        }
        cout << endl;
    }
    template <typename Seq>
    void prints(Seq &sq) {
        typename Seq::iterator it = sq.begin();
        while (it != sq.end())
            cout << *(it++) << "\t";
        cout << endl;
    }
    template <typename type>
    void prints(const deque<type> & sq) {
        for (unsigned long i = 0; i < sq.size(); i++)
            cout << sq[i] << "\t";
        cout << endl;
    }
    template<typename T, template<typename> class C>
    void printm(C<T>& c) {
        typename C<T>::iterator it = c.begin();
        while (it != c.end()) {
            prints(*it);
            it++;
        }
        cout << endl;
    }
    double ran2(long *idum) {
        long j;
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
    double ran4(bool t, long s) {
        double r = 0;
        static long seed_ = 1;
        if (t)
            r = ran2(&seed_);
        else
            seed_ = s;
        return r;
    }
    double ran4() {
        return ran4(true, 0);
    }
    void srand4(void) {
        long s = (long) time(NULL);
        ran4(false, s);
    }
    void srand5(long rank) {
        long s = (long) (rank);
        ran4(false, s);
    }
    long irand(long n) {
        return (int(ran4()*(n + 1)));
    }
    /*void srand_file(void) {
            ifstream in("time_seed.dat");
            long seed;
            if (!in.is_open())
                    seed=21111983;
            else
                    in>>seed;
            if (seed < 1 || seed>R2_IM2)
                    seed=1;
            srand5(seed);
            ofstream out("time_seed.dat");
            out<<seed+1<<endl;
    }*/
    void srand_file(void) {
        srand5(time(NULL));
    }
    template <typename Seq>
    double average_func(Seq &sq) {
        if (sq.empty())
            return 0;
        double av = 0;
        typename Seq::iterator it = sq.begin();
        while (it != sq.end())
            av += *(it++);
        av = av / sq.size();
        return av;
    }
    template <typename Seq>
    double variance_func(Seq &sq) {
        if (sq.empty())
            return 0;
        double av = 0;
        double var = 0;
        typename Seq::iterator it = sq.begin();
        while (it != sq.end()) {
            av += *(it);
            var += (*(it))*(*(it));
            it++;
        }
        av = av / sq.size();
        var = var / sq.size();
        var -= av*av;
        if (var < 1e-7)
            return 0;
        return var;
    }
    // this returns the average of the discrete probability function stored in Seq
    template <typename Seq>
    double average_pf(Seq &sq) {
        double av = 0;
        long h = 0;
        typename Seq::iterator it = sq.begin();
        while (it != sq.end()) {
            av += *(it) * h;
            it++;
            h++;
        }
        return av;
    }
    template <typename Seq>
    double variance_pf(Seq &sq) {
        double av = 0;
        double var = 0;
        long h = 0;
        typename Seq::iterator it = sq.begin();
        while (it != sq.end()) {
            av += *(it) * h;
            var += (*(it)) * h * h;
            it++;
            h++;
        }
        var -= av*av;
        if (var < 1e-7)
            return 0;
        return var;
    }
    double log_factorial(long num) {
        double log_result = 0;
        for (long i = 1; i <= num; i++)
            log_result += log(i);
        return (log_result);
    }

    double log_combination(long n, long k) {
        if (k == 0) return 0;
        if (n < k) return 0;
        if (n - k < k) k = n - k;
        double log_c = 0;
        for (long i = n - k + 1; i <= n; i++) log_c += log(i);
        for (long i = 1; i <= k; i++) log_c -= log(i);
        return log_c;
    }
    
    double binomial(long n, long x, double p) { //	returns the binomial distribution, n trials, x successes, p probability
        if (p == 0)
        {
            if (x == 0) return 1;
            else return 0;
        }
        if (p >= 1)
        {
            if (x == n) return 1;
            else return 0;
        }
        double log_b = 0;
        log_b += log_combination(n, x) + x * log(p)+(n - x) * log(1 - p);
        return (exp(log_b));
    }
    //to draw a number:
    /*
    deque <double> cumulative;
    binomial_cumulative(10, 0.5, cumulative);
    long nn=lower_bound(cumulative.begin(), cumulative.end(), ran4())-cumulative.begin();
     */
    long binomial_cumulative(long n, double p, deque<double> &cum) {
        cum.clear();
        double c = 0;
        for (long i = 0; i <= n; i++) {
            c += binomial(n, i, p);
            cum.push_back(c);
        }
        return 0;
    }
    // this function sets "cumulative" as the cumulative function of (1/x)^tau, with range= [min, n]
    //to draw a number:
    //long nn=lower_bound(cumulative.begin(), cumulative.end(), ran4())-cumulative.begin()+min_degree;
    long powerlaw(long n, long min, double tau, deque<double> &cumulative) {
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
    long distribution_from_cumulative(const deque<double> &cum, deque<double> &distr) { // cum is the cumulative, distr is set equal to the distribution
        distr.clear();
        double previous = 0;
        for (unsigned long i = 0; i < cum.size(); i++) {
            distr.push_back(cum[i] - previous);
            previous = cum[i];
        }
        return 0;
    }
    long cumulative_from_distribution(deque<double> &cum, const deque<double> &distr) { // cum is set equal to the cumulative, distr is the distribution
        cum.clear();
        double sum = 0;
        for (unsigned long i = 0; i < distr.size(); i++) {
            sum += distr[i];
            cum.push_back(sum);
        }
        return 0;
    }
    double poisson(long x, double mu) {
        return (exp(-mu + x * log(mu) - log_factorial(x)));
    }
    long shuffle_and_set(long *due, const long &dim) { // it sets due as a random sequence of integers from 0 to dim-1
        multimap <double, long> uno;
        for (long i = 0; i < dim; i++)
            uno.insert(make_pair(ran4(), i));
        multimap<double, long>::iterator it;
        long h = 0;
        for (it = uno.begin(); it != uno.end(); it++)
            due[h++] = it->second;
        return 0;
    }
    long shuffle_s(deque<long> &sq) {
        long siz = sq.size();
        if (siz == 0)
            return -1;
        for (unsigned long i = 0; i < sq.size(); i++) {
            long random_pos = irand(siz - 1);
            long random_card_ = sq[random_pos];
            sq[random_pos] = sq[siz - 1];
            sq[siz - 1] = random_card_;
            siz--;
        }
        return 0;
    }
    double compute_r(long x, long k, long kout, long m) {
        double r = 0;
        for (long i = x; i <= k; i++)
            r += binomial(k, i, double(kout) / double(m));
        return r;
    }
    long add_factors(deque<double> & num, deque<double> &den, long n, long k) {
        if (n < k)
            return -1;
        if (n - k < k)
            k = n - k;
        if (k == 0)
            return 0;
        for (long i = n - k + 1; i <= n; i++)
            num.push_back(double(i));
        for (long i = 1; i <= k; i++)
            den.push_back(double(i));
        return 0;
    }
    double compute_hypergeometric(long i, long k, long kout, long m) {
        if (i > k || i > kout || k > m || kout > m)
            return 0;
        double prod = 1;
        deque <double> num;
        deque <double> den;
        if (add_factors(num, den, kout, i) == -1)
            return 0;
        if (add_factors(num, den, m - kout, k - i) == -1)
            return 0;
        if (add_factors(den, num, m, k) == -1)
            return 0;
        sort(num.begin(), num.end());
        sort(den.begin(), den.end());
        //prints(den);
        for (unsigned long h = 0; h < den.size(); h++) if (den[h] <= 0) {
                cerr << "denominator has zero or less (in the hypergeometric)" << endl;
                return 0;
            }
        for (unsigned long h = 0; h < num.size(); h++) if (num[h] <= 0) {
                cerr << "numerator has zero or less (in the hypergeometric)" << endl;
                return 0;
            }
        //cout<<"sizes: "<<num.size()<<" "<<den.size()<<endl;
        for (unsigned long i = 0; i < num.size(); i++)
            prod = prod * num[i] / den[i];
        return prod;
    }
    //*
    double compute_self_links(long k, long n, long x) {
        if (2 * x > k)
            return 0;
        double prod = log_combination(n / 2, k - x) + log_combination(k - x, x) + (k - 2 * x) * log(2) - log_combination(n, k);
        return exp(prod);
    }
    //*/
    template <typename type>
    long log_histogram(deque<type> &c, ostream & out, long number_of_bins) { // c is the set od data, min is the lower bound, max is the upper one
        deque <type> d;
        for (unsigned long i = 0; i < c.size(); i++) if (c[i] > 0)
                d.push_back(c[i]);
        c.clear();
        c = d;
        double min = double(c[0]);
        double max = double(c[0]);
        for (unsigned long i = 0; i < c.size(); i++) {
            if (min>double(c[i]))
                min = double(c[i]);
            if (max<double(c[i]))
                max = double(c[i]);
        }
        deque <long> hist;
        deque <double> hist2;
        deque <double> bins;
        double step = log(min);
        if (max == min)
            max++;
        double bin = (log(max) - log(min)) / number_of_bins; // bin width
        while (step <= log(max) + 2 * bin) {
            bins.push_back(exp(step));
            hist.push_back(0);
            hist2.push_back(0);
            step += bin;
        }
        for (unsigned long i = 0; i < c.size(); i++) {
            long index = bins.size() - 1;
            for (unsigned long j = 0; j < bins.size() - 1; j++) if ((fabs(double(c[i]) - bins[j]) < 1e-7) || (double(c[i]) > bins[j] && double(c[i]) < bins[j + 1])) {
                    // this could be done in a more efficient way
                    index = j;
                    break;
                }
            //cout<<hist[index]<<" "<<index<<endl;
            hist[index]++;
            hist2[index] += double(c[i]);
        }
        for (unsigned long i = 0; i < hist.size() - 1; i++) {
            double h1 = bins[i];
            double h2 = bins[i + 1];
            double x = hist2[i] / hist[i];
            double y = double(hist[i]) / (c.size()*(h2 - h1));
            if (fabs(y) > 1e-10)
                out << x << "\t" << y << endl;
        }
        return 0;
    }
    template <typename type>
    long histogram(vector <type> &c, ostream & out, long number_of_bins, double b1, double b2) {
        // this should be OK
        // c is the set of data, b1 is the lower bound, b2 is the upper one (if they are equal, default limits are used)
        double min = double(c[0]);
        double max = double(c[0]);
        for (unsigned long i = 0; i < c.size(); i++) {
            if (min>double(c[i]))
                min = double(c[i]);
            if (max<double(c[i]))
                max = double(c[i]);
        }
        min -= 1e-6;
        max += 1e-6;
        if (b1 != b2) {
            min = b1;
            max = b2;
        }
        if (max == min)
            max += 1e-3;
        deque <long> hist;
        deque <double> hist2;
        double step = min;
        double bin = (max - min) / number_of_bins; // bin width
        while (step <= max + 2 * bin) {
            hist.push_back(0);
            hist2.push_back(0);
            step += bin;
        }
        for (unsigned long i = 0; i < c.size(); i++) {
            double data = double(c[i]);
            if (data > min && data <= max) {
                long index = int((data - min) / bin);
                hist[index]++;
                hist2[index] += double(c[i]);
            }
        }
        for (long i = 0; i < hist.size() - 1; i++) {
            double x = hist2[i] / hist[i];
            double y = double(hist[i]) / (c.size() * bin);
            if (fabs(y) > 1e-10)
                out << x << "\t" << y << endl;
        }
        return 0;
    }
    template <typename type>
    long not_norm_histogram_correlated(deque<type> &c, deque<type> &d, ostream & out, long number_of_bins, double b1, double b2) {
        // c is the x axis, d the y, b1 is the lower bound, b2 is the upper one (if they are equal, default limits are used)
        double min = double(c[0]);
        double max = double(c[0]);
        for (unsigned long i = 0; i < c.size(); i++) {
            if (min>double(c[i]))
                min = double(c[i]);
            if (max<double(c[i]))
                max = double(c[i]);
        }
        min -= 1e-6;
        max += 1e-6;
        if (b1 != b2) {
            min = b1;
            max = b2;
        }
        if (max == min)
            max += 1e-3;
        deque <long> hist; // frequency in the bin
        deque <double> hist_x; // x sum in the bin
        deque <double> hist_y; // y sum in the bin
        double step = min;
        double bin = (max - min) / number_of_bins; // bin width
        while (step <= max + 2 * bin) {
            hist.push_back(0);
            hist_x.push_back(0);
            hist_y.push_back(0);
            step += bin;
        }
        for (unsigned long i = 0; i < c.size(); i++) {
            double data = double(c[i]);
            if (data > min && data <= max) {
                long index = int((data - min) / bin);
                hist[index]++;
                hist_x[index] += double(c[i]);
                hist_y[index] += double(d[i]);
            }
        }
        for (long i = 0; i < hist.size() - 1; i++) {
            double x = hist_x[i] / hist[i];
            double y = hist_y[i] / hist[i];
            ;
            if (fabs(y) > 1e-10)
                out << x << "\t" << y << endl;
        }
        return 0;
    }
    template <typename type>
    long histogram(deque <type> &c, ostream & out, long number_of_bins, double b1, double b2) {
        // this should be OK
        // c is the set of data, b1 is the lower bound, b2 is the upper one (if they are equal, default limits are used)
        double min = double(c[0]);
        double max = double(c[0]);
        for (unsigned long i = 0; i < c.size(); i++) {
            if (min>double(c[i]))
                min = double(c[i]);
            if (max<double(c[i]))
                max = double(c[i]);
        }
        min -= 1e-6;
        max += 1e-6;
        if (b1 != b2) {
            min = b1;
            max = b2;
        }
        if (max == min)
            max += 1e-3;
        deque <long> hist;
        deque <double> hist2;
        double step = min;
        double bin = (max - min) / number_of_bins; // bin width
        while (step <= max + 2 * bin) {
            hist.push_back(0);
            hist2.push_back(0);
            step += bin;
        }
        for (unsigned long i = 0; i < c.size(); i++) {
            double data = double(c[i]);
            if (data > min && data <= max) {
                long index = int((data - min) / bin);
                hist[index]++;
                hist2[index] += double(c[i]);
            }
        }
        for (long i = 0; i < hist.size() - 1; i++) {
            double x = hist2[i] / hist[i];
            double y = double(hist[i]) / (c.size() * bin);
            if (fabs(y) > 1e-10)
                out << x << "\t" << y << endl;
        }
        return 0;
    }
    template <typename type>
    long not_norm_histogram(vector<type> &c, ostream & out, long number_of_bins, double b1, double b2) {
        // this should be OK
        // c is the set of data, b1 is the lower bound, b2 is the upper one (if they are equal, default limits are used)
        double min = double(c[0]);
        double max = double(c[0]);
        for (unsigned long i = 0; i < c.size(); i++) {
            if (min>double(c[i]))
                min = double(c[i]);
            if (max<double(c[i]))
                max = double(c[i]);
        }
        min -= 1e-6;
        max += 1e-6;
        if (b1 != b2) {
            min = b1;
            max = b2;
        }
        if (max == min)
            max += 1e-3;
        deque <long> hist;
        deque <double> hist2;
        double step = min;
        double bin = (max - min) / number_of_bins; // bin width
        while (step <= max + 2 * bin) {
            hist.push_back(0);
            hist2.push_back(0);
            step += bin;
        }
        for (unsigned long i = 0; i < c.size(); i++) {
            double data = double(c[i]);
            if (data > min && data <= max) {
                long index = int((data - min) / bin);
                hist[index]++;
                hist2[index] += double(c[i]);
            }
        }
        for (long i = 0; i < hist.size() - 1; i++) {
            double x = hist2[i] / hist[i];
            double y = double(hist[i]) / (c.size());
            if (fabs(y) > 1e-10)
                out << x << "\t" << y << endl;
        }
        return 0;
    }
    template <typename type>
    long not_norm_histogram(deque<type> &c, ostream & out, long number_of_bins, double b1, double b2) {
        // this should be OK
        // c is the set of data, b1 is the lower bound, b2 is the upper one (if they are equal, default limits are used)
        double min = double(c[0]);
        double max = double(c[0]);
        for (unsigned long i = 0; i < c.size(); i++) {
            if (min>double(c[i]))
                min = double(c[i]);
            if (max<double(c[i]))
                max = double(c[i]);
        }
        min -= 1e-6;
        max += 1e-6;
        if (b1 != b2) {
            min = b1;
            max = b2;
        }
        if (max == min)
            max += 1e-3;
        deque <long> hist;
        deque <double> hist2;
        double step = min;
        double bin = (max - min) / number_of_bins; // bin width
        while (step <= max + 2 * bin) {
            hist.push_back(0);
            hist2.push_back(0);
            step += bin;
        }
        for (unsigned long i = 0; i < c.size(); i++) {
            double data = double(c[i]);
            if (data > min && data <= max) {
                long index = int((data - min) / bin);
                hist[index]++;
                hist2[index] += double(c[i]);
            }
        }
        for (unsigned long i = 0; i < hist.size() - 1; i++) {
            double x = hist2[i] / hist[i];
            double y = double(hist[i]) / (c.size());
            if (fabs(y) > 1e-10)
                out << x << "\t" << y << endl;
        }
        return 0;
    }
    long int_histogram(vector <long> &c, ostream & out) {
        map<long, double> hist;
        double freq = 1 / double(c.size());
        for (unsigned long i = 0; i < c.size(); i++) {
            map<long, double>::iterator itf = hist.find(c[i]);
            if (itf == hist.end())
                hist.insert(make_pair(c[i], 1.));
            else
                itf->second++;
        }
        for (map<long, double>::iterator it = hist.begin(); it != hist.end(); it++)
            it->second = it->second * freq;
        prints(hist, out);
        return 1;
    }
    long int_histogram(deque <long> &c, ostream & out) {
        map<long, double> hist;
        double freq = 1 / double(c.size());
        for (unsigned long i = 0; i < c.size(); i++) {
            map<long, double>::iterator itf = hist.find(c[i]);
            if (itf == hist.end())
                hist.insert(make_pair(c[i], 1.));
            else
                itf->second++;
        }
        for (map<long, double>::iterator it = hist.begin(); it != hist.end(); it++)
            it->second = it->second * freq;
        prints(hist, out);
        return 1;
    }
    bool cast_string_to_double(string &b, double &h) {
        // set h= the number written in b[];
        // return false if there is an error
        h = 0;
        if (b.size() == 0)
            return false;
        long sign = 1;
        if (b[0] == '-') {
            b[0] = '0';
            sign = -1;
        }
        unsigned long digits_before = 0;
        for (unsigned long i = 0; i < b.size(); i++)
            if (b[i] != '.')
                digits_before++;
            else
                break;
        unsigned long j = 0;
        while (j != digits_before) {
            long number = (int(b[j]) - 48);
            h += number * pow(10, digits_before - j - 1);
            if (number < 0 || number > 9)
                return false;
            j++;
        }
        j = digits_before + 1;
        while (j < b.size()) {
            long number = (int(b[j]) - 48);
            h += number * pow(10, digits_before - j);
            if (number < 0 || number > 9)
                return false;
            j++;
        }
        h = sign*h;
        return true;
    }
    long cast_int(double u) {
        long a = int(u);
        if (u - a > 0.5)
            a++;
        return a;
    }
    long cast_string_to_char(string &file_name, char *b) {
        for (unsigned long i = 0; i < file_name.size(); i++)
            b[i] = file_name[i];
        b[file_name.size()] = '\0';
        return 0;
    }
    class Parameters {
    public:
        Parameters();
        ~Parameters() {
        };
        long num_nodes;
        double average_k;
        long max_degree;
        double tau;
        double tau2;
        double mixing_parameter;
        long overlapping_nodes;
        long overlap_membership;
        long nmin;
        long nmax;
        bool fixed_range;
        bool excess;
        bool defect;
        bool randomf;
        bool set(string &, string &);
        void set_random();
        bool arrange();
        deque<string> command_flags;
    };
    Parameters::Parameters() {
        num_nodes = unlikely;
        average_k = unlikely;
        max_degree = unlikely;
        tau = 2;
        tau2 = 1;
        mixing_parameter = unlikely;
        overlapping_nodes = 0;
        overlap_membership = 0;
        nmin = unlikely;
        nmax = unlikely;
        randomf = false;
        fixed_range = false;
        excess = false;
        defect = false;
        command_flags.push_back("-N"); //0
        command_flags.push_back("-k"); //1
        command_flags.push_back("-maxk"); //2
        command_flags.push_back("-mu"); //3
        command_flags.push_back("-t1"); //4
        command_flags.push_back("-t2"); //5
        command_flags.push_back("-minc"); //6
        command_flags.push_back("-maxc"); //7
        command_flags.push_back("-on"); //8
        command_flags.push_back("-om"); //9
    };
    void Parameters::set_random() {
        cout << "this is a random network" << endl;
        mixing_parameter = 0;
        overlapping_nodes = 0;
        overlap_membership = 0;
        nmax = num_nodes;
        nmin = num_nodes;
        fixed_range = true;
        excess = false;
        defect = false;
    }
    bool Parameters::arrange() {
        if (randomf)
            set_random();
        if (num_nodes == unlikely) {
            cerr << "\n***********************\nERROR:\t number of nodes unspecified" << endl;
            return false;
        }
        if (average_k == unlikely) {
            cerr << "\n***********************\nERROR:\t average degree unspecified" << endl;
            return false;
        }
        if (max_degree == unlikely) {
            cerr << "\n***********************\nERROR:\t maximum degree unspecified" << endl;
            return false;
        }
        if (mixing_parameter == unlikely) {
            cerr << "\n***********************\nERROR:\t mixing parameter unspecified" << endl;
            return false;
        }
        if (overlapping_nodes < 0 || overlap_membership < 0) {
            cerr << "\n***********************\nERROR:\tsome positive parameters are negative" << endl;
            return -1;
        }
        if (num_nodes <= 0 || average_k <= 0 || max_degree <= 0 || mixing_parameter < 0 || (nmax <= 0 && nmax != unlikely) || (nmin <= 0 && nmin != unlikely)) {
            cerr << "\n***********************\nERROR:\tsome positive parameters are negative" << endl;
            return -1;
        }
        if (mixing_parameter > 1) {
            cerr << "\n***********************\nERROR:\tmixing parameter > 1 (must be between 0 and 1)" << endl;
            return -1;
        }
        if (nmax != unlikely && nmin != unlikely)
            fixed_range = true;
        else
            fixed_range = false;
        if (excess && defect) {
            cerr << "\n***********************\nERROR:\tboth options -inf and -sup cannot be used at the same time" << endl;
            return false;
        }
        cout << "\n**************************************************************" << endl;
        cout << "number of nodes:\t" << num_nodes << endl;
        cout << "average degree:\t" << average_k << endl;
        cout << "maximum degree:\t" << max_degree << endl;
        cout << "exponent for the degree distribution:\t" << tau << endl;
        cout << "exponent for the community size distribution:\t" << tau2 << endl;
        cout << "mixing parameter:\t" << mixing_parameter << endl;
        cout << "number of overlapping nodes:\t" << overlapping_nodes << endl;
        cout << "number of memberships of the overlapping nodes:\t" << overlap_membership << endl;
        if (fixed_range) {
            cout << "community size range set equal to [" << nmin << " , " << nmax << "]" << endl;
            if (nmin > nmax) {
                cerr << "\n***********************\nERROR: INVERTED COMMUNITY SIZE BOUNDS" << endl;
                return false;
            }
            if (nmax > num_nodes) {
                cerr << "\n***********************\nERROR: maxc BIGGER THAN THE NUMBER OF NODES" << endl;
                return false;
            }
        }
        cout << "**************************************************************" << endl << endl;
        return true;
    }
    bool Parameters::set(string & flag, string & num) {
        // false is something goes wrong
        cout << "setting... " << flag << " " << num << endl;
        double err;
        if (!cast_string_to_double(num, err)) {
            cerr << "\n***********************\nERROR while reading parameters" << endl;
            return false;
        }
        if (flag == command_flags[0]) {
            if (fabs(err - long (err)) > 1e-8) {
                cerr << "\n***********************\nERROR: number of nodes must be an integer" << endl;
                return false;
            }
            num_nodes = cast_int(err);
        } else if (flag == command_flags[1]) {
            average_k = err;
        } else if (flag == command_flags[2]) {
            max_degree = cast_int(err);
        } else if (flag == command_flags[3]) {
            mixing_parameter = err;
        } else if (flag == command_flags[4]) {
            tau = err;
        } else if (flag == command_flags[5]) {
            tau2 = err;
        }
        else if (flag == command_flags[6]) {
            if (fabs(err - long (err)) > 1e-8) {
                cerr << "\n***********************\nERROR: the minumum community size must be an integer" << endl;
                return false;
            }
            nmin = cast_int(err);
        }
        else if (flag == command_flags[7]) {
            if (fabs(err - long (err)) > 1e-8) {
                cerr << "\n***********************\nERROR: the maximum community size must be an integer" << endl;
                return false;
            }
            nmax = cast_int(err);
        } else if (flag == command_flags[8]) {
            if (fabs(err - long (err)) > 1e-8) {
                cerr << "\n***********************\nERROR: the number of overlapping nodes must be an integer" << endl;
                return false;
            }
            overlapping_nodes = cast_int(err);
        } else if (flag == command_flags[9]) {
            if (fabs(err - long (err)) > 1e-8) {
                cerr << "\n***********************\nERROR: the number of membership of the overlapping nodes must be an integer" << endl;
                return false;
            }
            overlap_membership = cast_int(err);
        } else {
            cerr << "\n***********************\nERROR while reading parameters: " << flag << " is an unknown option" << endl;
            return false;
        }
        return true;
    }
    void statement() {
        cout << "\nTo run the program type \n./benchmark [FLAG] [P]" << endl;
        cout << "\n----------------------\n" << endl;
        cout << "To set the parameters, type:" << endl << endl;
        cout << "-N\t\t[number of nodes]" << endl;
        cout << "-k\t\t[average degree]" << endl;
        cout << "-maxk\t\t[maximum degree]" << endl;
        cout << "-mu\t\t[mixing parameter]" << endl;
        cout << "-t1\t\t[minus exponent for the degree sequence]" << endl;
        cout << "-t2\t\t[minus exponent for the community size distribution]" << endl;
        cout << "-minc\t\t[minimum for the community sizes]" << endl;
        cout << "-maxc\t\t[maximum for the community sizes]" << endl;
        cout << "-on\t\t[number of overlapping nodes]" << endl;
        cout << "-om\t\t[number of memberships of the overlapping nodes]" << endl;
        cout << "----------------------\n" << endl;
        cout << "It is also possible to set the parameters writing flags and relative numbers in a file. To specify the file, use the option:" << endl;
        cout << "-f\t[filename]" << endl;
        cout << "You can set the parameters both writing some of them in the file, and using flags from the command line for others." << endl << endl;
        cout << "-N, -k, -maxk, -mu have to be specified. For the others, the program can use default values:" << endl;
        cout << "t1=2, t2=1, on=0, om=0, minc and maxc will be chosen close to the degree sequence extremes." << endl;
        cout << "If you set a parameter twice, the latter one will be taken." << endl;
        cout << "\n-------------------- Other options ---------------------------\n" << endl;
        cout << "To have a random network use:" << endl;
        cout << "-rand" << endl;
        cout << "Using this option will set mu=0, and minc=maxc=N, i.e. there will be one only community." << endl;
        cout << "Use option -sup (-inf) if you want to produce a benchmark whose distribution of the ratio of external degree/total degree ";
        cout << "is superiorly (inferiorly) bounded by the mixing parameter." << endl;
        cout << "\n-------------------- Examples ---------------------------\n" << endl;
        cout << "Example1:" << endl;
        cout << "./benchmark -N 1000 -k 15 -maxk 50 -mu 0.1 -minc 20 -maxc 50" << endl;
        cout << "Example2:" << endl;
        cout << "./benchmark -f flags.dat -t1 3" << endl;
        cout << "\n-------------------- Other info ---------------------------\n" << endl;
        cout << "Read file ReadMe.txt for more info." << endl << endl;
    }
    bool set_from_file(string & file_name, Parameters & par1) {
        long h = file_name.size();
        char b[h + 1];
        cast_string_to_char(file_name, b);
        ifstream in(b);
        if (!in.is_open()) {
            cerr << "File " << file_name << " not found. Where is it?" << endl;
            return false;
        }
        string temp;
        while (in >> temp) { // input file name
            if (temp == "-rand")
                par1.randomf = true;
            else if (temp == "-sup")
                par1.excess = true;
            else if (temp == "-inf")
                par1.defect = true;
            else {
                string temp2;
                in >> temp2;
                if (temp2.size() > 0) {
                    if (temp == "-f" && temp2 != file_name) {
                        if (set_from_file(temp2, par1) == false)
                            return false;
                    }
                    if (temp != "-f") {
                        if (par1.set(temp, temp2) == false)
                            return false;
                    }
                }
                else {
                    cerr << "\n***********************\nERROR while reading parameters" << endl;
                    return false;
                }
            }
        }
        return true;
    }
    bool set_parameters(long argc, char * argv[], Parameters & par1) {
        long argct = 0;
        string temp;
        if (argc <= 1) { // if no arguments, return statement about program usage.
            statement();
            return false;
        }
        while (++argct < argc) { // input file name
            temp = argv[argct];
            if (temp == "-rand")
                par1.randomf = true;
            else if (temp == "-sup")
                par1.excess = true;
            else if (temp == "-inf")
                par1.defect = true;
            else {
                argct++;
                string temp2;
                if (argct < argc) {
                    temp2 = argv[argct];
                    if (temp == "-f") {
                        if (set_from_file(temp2, par1) == false)
                            return false;
                    }
                    if (temp != "-f") {
                        if (par1.set(temp, temp2) == false)
                            return false;
                    }
                }
                else {
                    cerr << "\n***********************\nERROR while reading parameters" << endl;
                    return false;
                }
            }
        }
        if (par1.arrange() == false)
            return false;
        return true;
    }
    // it computes the sum of a deque<long>
    long deque_int_sum(const deque<long> & a) {
        long s = 0;
        for (unsigned long i = 0; i < a.size(); i++)
            s += a[i];
        return s;
    }
    // it computes the integral of a power law
    double integral(double a, double b) {
        if (fabs(a + 1.) > 1e-10)
            return (1. / (a + 1.) * pow(b, a + 1.));
        else
            return (log(b));
    }
    // it returns the average degree of a power law
    double average_degree(const double &dmax, const double &dmin, const double &gamma) {
        return (1. / (integral(gamma, dmax) - integral(gamma, dmin)))*(integral(gamma + 1, dmax) - integral(gamma + 1, dmin));
    }
    //bisection method to find the inferior limit, in order to have the expected average degree
    double solve_dmin(const double& dmax, const double &dmed, const double &gamma) {
        double dmin_l = 1;
        double dmin_r = dmax;
        double average_k1 = average_degree(dmin_r, dmin_l, gamma);
        double average_k2 = dmin_r;
        if ((average_k1 - dmed > 0) || (average_k2 - dmed < 0)) {
            cerr << "\n***********************\nERROR: the average degree is out of range:";
            if (average_k1 - dmed > 0) {
                cerr << "\nyou should increase the average degree (bigger than " << average_k1 << ")" << endl;
                cerr << "(or decrease the maximum degree...)" << endl;
            }
            if (average_k2 - dmed < 0) {
                cerr << "\nyou should decrease the average degree (smaller than " << average_k2 << ")" << endl;
                cerr << "(or increase the maximum degree...)" << endl;
            }
            return -1;
        }
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
    // it computes the correct (i.e. discrete) average of a power law
    double integer_average(long n, long min, double tau) {
        double a = 0;
        for (double h = min; h < n + 1; h++)
            a += pow((1. / h), tau);
        double pf = 0;
        for (double i = min; i < n + 1; i++)
            pf += 1 / a * pow((1. / (i)), tau) * i;
        return pf;
    }
    // this function changes the community sizes merging the smallest communities
    long change_community_size(deque<long> &seq) {
        if (seq.size() <= 2)
            return -1;
        unsigned long min1 = 0;
        unsigned long min2 = 0;
        for (unsigned long i = 0; i < seq.size(); i++)
            if (seq[i] <= seq[min1])
                min1 = i;
        if (min1 == 0)
            min2 = 1;
        for (unsigned long i = 0; i < seq.size(); i++)
            if (seq[i] <= seq[min2] && seq[i] > seq[min1])
                min2 = i;
        seq[min1] += seq[min2];
        long c = seq[0];
        seq[0] = seq[min2];
        seq[min2] = c;
        seq.pop_front();
        return 0;
    }
    long build_bipartite_network(deque<deque<long> > & member_matrix, const deque<long> & member_numbers, const deque<long> &num_seq) {
        // this function builds a bipartite network with num_seq and member_numbers which are the degree sequences. in member matrix links of the communities are stored
        // this means member_matrix has num_seq.size() rows and each row has num_seq[i] elements
        deque<set<long> > en_in; // this is the Ein of the subgraph
        deque<set<long> > en_out; // this is the Eout of the subgraph
        {
            set<long> first;
            for (unsigned long i = 0; i < member_numbers.size(); i++) {
                en_in.push_back(first);
            }
        }
        {
            set<long> first;
            for (unsigned long i = 0; i < num_seq.size(); i++) {
                en_out.push_back(first);
            }
        }
        multimap <long, long> degree_node_out;
        deque<pair<long, long> > degree_node_in;
        for (unsigned long i = 0; i < num_seq.size(); i++)
            degree_node_out.insert(make_pair(num_seq[i], i));
        for (unsigned long i = 0; i < member_numbers.size(); i++)
            degree_node_in.push_back(make_pair(member_numbers[i], i));
        sort(degree_node_in.begin(), degree_node_in.end());
        deque<pair<long, long> >::iterator itlast = degree_node_in.end();
        /*
        for (long i=0; i<degree_node_in.size(); i++)
                cout<<degree_node_in[i].first<<" "<<degree_node_in[i].second<<endl;
         */
        while (itlast != degree_node_in.begin()) {
            itlast--;
            multimap <long, long>::iterator itit = degree_node_out.end();
            deque <multimap<long, long>::iterator> erasenda;
            for (long i = 0; i < itlast->first; i++) {
                if (itit != degree_node_out.begin()) {
                    itit--;
                    en_in[itlast->second].insert(itit->second);
                    en_out[itit->second].insert(itlast->second);
                    erasenda.push_back(itit);
                }
                else
                    return -1;
            }
            //cout<<"degree node out before"<<endl;
            //prints(degree_node_out);
            for (unsigned long i = 0; i < erasenda.size(); i++) {
                if (erasenda[i]->first > 1)
                    degree_node_out.insert(make_pair(erasenda[i]->first - 1, erasenda[i]->second));
                degree_node_out.erase(erasenda[i]);
            }
            //cout<<"degree node out after"<<endl;
            //prints(degree_node_out);
        }
        // this is to randomize the subgraph -------------------------------------------------------------------
        for (unsigned long node_a = 0; node_a < num_seq.size(); node_a++) for (unsigned long krm = 0; krm < en_out[node_a].size(); krm++) {
                long random_mate = irand(member_numbers.size() - 1);
                if (en_out[node_a].find(random_mate) == en_out[node_a].end()) {
                    deque <long> external_nodes;
                    for (set<long>::iterator it_est = en_out[node_a].begin(); it_est != en_out[node_a].end(); it_est++)
                        external_nodes.push_back(*it_est);
                    long old_node = external_nodes[irand(external_nodes.size() - 1)];
                    deque <long> not_common;
                    for (set<long>::iterator it_est = en_in[random_mate].begin(); it_est != en_in[random_mate].end(); it_est++)
                        if (en_in[old_node].find(*it_est) == en_in[old_node].end())
                            not_common.push_back(*it_est);
                    if (not_common.empty())
                        break;
                    long node_h = not_common[irand(not_common.size() - 1)];
                    en_out[node_a].insert(random_mate);
                    en_out[node_a].erase(old_node);
                    en_in[old_node].insert(node_h);
                    en_in[old_node].erase(node_a);
                    en_in[random_mate].insert(node_a);
                    en_in[random_mate].erase(node_h);
                    en_out[node_h].erase(random_mate);
                    en_out[node_h].insert(old_node);
                }
            }
        member_matrix.clear();
        deque <long> first;
        for (unsigned long i = 0; i < en_out.size(); i++) {
            member_matrix.push_back(first);
            for (set<long>::iterator its = en_out[i].begin(); its != en_out[i].end(); its++)
                member_matrix[i].push_back(*its);
        }
        return 0;
    }
    long internal_degree_and_membership(double mixing_parameter, long overlapping_nodes, long max_mem_num, long num_nodes, deque<deque<long> > & member_matrix,
            bool excess, bool defect, deque<long> & degree_seq, deque<long> &num_seq, deque<long> &internal_degree_seq, bool fixed_range, long nmin, long nmax, double tau2) {
        if (num_nodes < overlapping_nodes) {
            cerr << "\n***********************\nERROR: there are more overlapping nodes than nodes in the whole network! Please, decrease the former ones or increase the latter ones" << endl;
            return -1;
        }
        //
        member_matrix.clear();
        internal_degree_seq.clear();
        deque<double> cumulative;
        // it assigns the internal degree to each node -------------------------------------------------------------------------
        long max_degree_actual = 0; // maximum internal degree
        for (unsigned long i = 0; i < degree_seq.size(); i++) {
            double interno = (1 - mixing_parameter) * degree_seq[i];
            long int_interno = int(interno);
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
        }
        // it assigns the community size sequence -----------------------------------------------------------------------------
        powerlaw(nmax, nmin, tau2, cumulative);
        if (num_seq.empty()) {
            long _num_ = 0;
            if (!fixed_range && (max_degree_actual + 1) > nmin) {
                _num_ = max_degree_actual + 1; // this helps the assignment of the memberships (it assures that at least one module is big enough to host each node)
                num_seq.push_back(max_degree_actual + 1);
            }
            while (true) {
                long nn = lower_bound(cumulative.begin(), cumulative.end(), ran4()) - cumulative.begin() + nmin;
                if (nn + _num_ <= num_nodes + overlapping_nodes * (max_mem_num - 1)) {
                    num_seq.push_back(nn);
                    _num_ += nn;
                } else
                    break;
            }
            num_seq[min_element(num_seq.begin(), num_seq.end()) - num_seq.begin()] += num_nodes + overlapping_nodes * (max_mem_num - 1) - _num_;
        }
        //cout<<"num_seq"<<endl;
        //prints(num_seq);
        //cout<<"\n----------------------------------------------------------"<<endl;
        /*
        cout<<"community sizes"<<endl;
        for (long i=0; i<num_seq.size(); i++)
                cout<<num_seq[i]<<" ";
        cout<<endl<<endl;
        //*/
        /*
        deque <long> first;
        for (long i=0; i<ncom; i++)
                member_matrix.push_back(first);
        // it puts the overlapping_nodes inside
        cout<<ncom<<endl;
        for (long i=degree_seq.size() - overlapping_nodes; i<degree_seq.size(); i++) {
                cout<<i<<endl;
                set<long> members;
                long hh=0;
                while(members.size()<max_mem_num) {
                        long random_module=irand(ncom-1);
                        if(member_matrix[random_module].size()!=num_seq[random_module])
                                members.insert(random_module);
                        hh++;
                        if(hh>3*num_nodes) {
                                cerr<<"it seems that the overlapping nodes need more communities that those I provided. Please increase the number of communities or decrease the number of overlapping nodes"<<endl;
                                return -1;
                        }
                }
                for (set<long>::iterator its=members.begin(); its!=members.end(); its++)
                        member_matrix[*its].push_back(i);
        }
        // it decides the memberships for the not overlapping nodes
        long moment_module=0;
        for (long i=0; i<num_nodes - overlapping_nodes; i++) {
                while(member_matrix[moment_module].size()==num_seq[moment_module])
                         moment_module++;
                member_matrix[moment_module].push_back(i);
        }
         */
        // I have to assign the degree to the nodes
        deque<long> member_numbers;
        for (long i = 0; i < overlapping_nodes; i++)
            member_numbers.push_back(max_mem_num);
        for (unsigned long i = overlapping_nodes; i < degree_seq.size(); i++)
            member_numbers.push_back(1);
        //prints(member_numbers);
        //prints(num_seq);
        if (build_bipartite_network(member_matrix, member_numbers, num_seq) == -1) {
            cerr << "it seems that the overlapping nodes need more communities that those I provided. Please increase the number of communities or decrease the number of overlapping nodes" << endl;
            return -1;
        }
        //printm(member_matrix);
        //cout<<"degree_seq"<<endl;
        //prints(degree_seq);
        //cout<<"internal_degree_seq"<<endl;
        //prints(internal_degree_seq);
        deque<long> available;
        for (long i = 0; i < num_nodes; i++)
            available.push_back(0);
        for (unsigned long i = 0; i < member_matrix.size(); i++) {
            for (unsigned long j = 0; j < member_matrix[i].size(); j++)
                available[member_matrix[i][j]] += member_matrix[i].size() - 1;
        }
        //cout<<"available"<<endl;
        //prints(available);
        deque<long> available_nodes;
        for (long i = 0; i < num_nodes; i++)
            available_nodes.push_back(i);
        deque<long> map_nodes; // in the position i there is the new name of the node i
        for (long i = 0; i < num_nodes; i++)
            map_nodes.push_back(0);
        for (long i = degree_seq.size() - 1; i >= 0; i--) {
            long try_this = irand(available_nodes.size() - 1);
            long kr = 0;
            while (internal_degree_seq[i] > available[available_nodes[try_this]]) {
                kr++;
                try_this = irand(available_nodes.size() - 1);
                if (kr == 3 * num_nodes) {
                    if (change_community_size(num_seq) == -1) {
                        cerr << "\n***********************\nERROR: this program needs more than one community to work fine" << endl;
                        return -1;
                    }
                    cout << "it took too long to decide the memberships; I will try to change the community sizes" << endl;
                    cout << "new community sizes" << endl;
                    for (unsigned long i = 0; i < num_seq.size(); i++)
                        cout << num_seq[i] << " ";
                    cout << endl << endl;
                    return (internal_degree_and_membership(mixing_parameter, overlapping_nodes, max_mem_num, num_nodes, member_matrix, excess, defect, degree_seq, num_seq, internal_degree_seq, fixed_range, nmin, nmax, tau2));
                }
            }
            map_nodes[available_nodes[try_this]] = i;
            available_nodes[try_this] = available_nodes[available_nodes.size() - 1];
            available_nodes.pop_back();
        }
        for (unsigned long i = 0; i < member_matrix.size(); i++) {
            for (unsigned long j = 0; j < member_matrix[i].size(); j++)
                member_matrix[i][j] = map_nodes[member_matrix[i][j]];
        }
        for (unsigned long i = 0; i < member_matrix.size(); i++)
            sort(member_matrix[i].begin(), member_matrix[i].end());
        return 0;
    }
    long compute_internal_degree_per_node(long d, long m, deque<long> & a) {
        // d is the internal degree
        // m is the number of memebership
        a.clear();
        long d_i = d / m;
        for (long i = 0; i < m; i++)
            a.push_back(d_i);
        for (long i = 0; i < d % m; i++)
            a[i]++;
        return 0;
    }
    /*
    long check_link_list(const deque<deque<long> > & link_list, const deque<long> & degree_seq) {
            for (long i=0; i<link_list.size(); i++) {
                    long s=0;
                    for (long j=0; j<link_list[i].size(); j++)
                            s+=link_list[i][j];
                    if(s!=degree_seq[i]) {
                            long ok;
                            cerr<<"wrong link list"<<endl;
                            cin>>ok;
                    }
            }
    }
     */
    long build_subgraph(deque<set<long> > & E, const deque<long> & nodes, const deque<long> & degrees) {
        /*
        cout<<"nodes"<<endl;
        prints(nodes);
        cout<<"degrees"<<endl;
        prints(degrees);
         */
        if (degrees.size() < 3) {
            cerr << "it seems that some communities should have only 2 nodes! This does not make much sense (in my opinion) Please change some parameters!" << endl;
            return -1;
        }
        // this function is to build a network with the labels stored in nodes and the degree seq in degrees (correspondence is based on the vectorial index)
        // the only complication is that you don't want the nodes to have neighbors they already have
        // labels will be placed in the end
        deque<set<long> > en; // this is the E of the subgraph
        {
            set<long> first;
            for (unsigned long i = 0; i < nodes.size(); i++)
                en.push_back(first);
        }
        multimap <long, long> degree_node;
        for (unsigned long i = 0; i < degrees.size(); i++)
            degree_node.insert(degree_node.end(), make_pair(degrees[i], i));
        long var = 0;
        while (degree_node.size() > 0) {
            multimap<long, long>::iterator itlast = degree_node.end();
            itlast--;
            multimap <long, long>::iterator itit = itlast;
            deque <multimap<long, long>::iterator> erasenda;
            long inserted = 0;
            for (long i = 0; i < itlast->first; i++) {
                if (itit != degree_node.begin()) {
                    itit--;
                    en[itlast->second].insert(itit->second);
                    en[itit->second].insert(itlast->second);
                    inserted++;
                    erasenda.push_back(itit);
                }
                else
                    break;
            }
            for (unsigned long i = 0; i < erasenda.size(); i++) {
                if (erasenda[i]->first > 1)
                    degree_node.insert(make_pair(erasenda[i]->first - 1, erasenda[i]->second));
                degree_node.erase(erasenda[i]);
            }
            var += itlast->first - inserted;
            degree_node.erase(itlast);
        }
        // this is to randomize the subgraph -------------------------------------------------------------------
        for (unsigned long node_a = 0; node_a < degrees.size(); node_a++) for (unsigned long krm = 0; krm < en[node_a].size(); krm++) {
                unsigned long random_mate = irand(degrees.size() - 1);
                while (random_mate == node_a)
                    random_mate = irand(degrees.size() - 1);
                if (en[node_a].insert(random_mate).second) {
                    deque <long> out_nodes;
                    for (set<long>::iterator it_est = en[node_a].begin(); it_est != en[node_a].end(); it_est++) if ((*it_est) != random_mate)
                            out_nodes.push_back(*it_est);
                    long old_node = out_nodes[irand(out_nodes.size() - 1)];
                    en[node_a].erase(old_node);
                    en[random_mate].insert(node_a);
                    en[old_node].erase(node_a);
                    deque <long> not_common;
                    for (set<long>::iterator it_est = en[random_mate].begin(); it_est != en[random_mate].end(); it_est++)
                        if ((old_node != (*it_est)) && (en[old_node].find(*it_est) == en[old_node].end()))
                            not_common.push_back(*it_est);
                    long node_h = not_common[irand(not_common.size() - 1)];
                    en[random_mate].erase(node_h);
                    en[node_h].erase(random_mate);
                    en[node_h].insert(old_node);
                    en[old_node].insert(node_h);
                }
            }
        // now I try to insert the new links into the already done network. If some multiple links come out, I try to rewire them
        deque < pair<long, long> > multiple_edge;
        for (unsigned long i = 0; i < en.size(); i++) {
            for (set<long>::iterator its = en[i].begin(); its != en[i].end(); its++) if (i<static_cast<unsigned long>(*its)) {
                    bool already = !(E[nodes[i]].insert(nodes[*its]).second); // true is the insertion didn't take place
                    if (already)
                        multiple_edge.push_back(make_pair(nodes[i], nodes[*its]));
                    else
                        E[nodes[*its]].insert(nodes[i]);
                }
        }
        //cout<<"multiples "<<multiple_edge.size()<<endl;
        for (unsigned long i = 0; i < multiple_edge.size(); i++) {
            long &a = multiple_edge[i].first;
            long &b = multiple_edge[i].second;
            // now, I'll try to rewire this multiple link among the nodes stored in nodes.
            unsigned long stopper_ml = 0;
            while (true) {
                stopper_ml++;
                long random_mate = nodes[irand(degrees.size() - 1)];
                while (random_mate == a || random_mate == b)
                    random_mate = nodes[irand(degrees.size() - 1)];
                if (E[a].find(random_mate) == E[a].end()) {
                    deque <long> not_common;
                    for (set<long>::iterator it_est = E[random_mate].begin(); it_est != E[random_mate].end(); it_est++)
                        if ((b != (*it_est)) && (E[b].find(*it_est) == E[b].end()) && (binary_search(nodes.begin(), nodes.end(), *it_est)))
                            not_common.push_back(*it_est);
                    if (not_common.size() > 0) {
                        long node_h = not_common[irand(not_common.size() - 1)];
                        E[random_mate].insert(a);
                        E[random_mate].erase(node_h);
                        E[node_h].erase(random_mate);
                        E[node_h].insert(b);
                        E[b].insert(node_h);
                        E[a].insert(random_mate);
                        break;
                    }
                }
                if (stopper_ml == 2 * E.size()) {
                    cout << "sorry, I need to change the degree distribution a little bit (one less link)" << endl;
                    break;
                }
            }
        }
        return 0;
    }
    long build_subgraphs(deque<set<long> > & E, const deque<deque<long> > & member_matrix, deque<deque<long> > & member_list,
            deque<deque<long> > & link_list, const deque<long> & internal_degree_seq, const deque<long> & degree_seq, const bool excess, const bool defect) {
        E.clear();
        member_list.clear();
        link_list.clear();
        long num_nodes = degree_seq.size();
        //printm(member_matrix);
        {
            deque<long> first;
            for (long i = 0; i < num_nodes; i++)
                member_list.push_back(first);
        }
        for (unsigned long i = 0; i < member_matrix.size(); i++)
            for (unsigned long j = 0; j < member_matrix[i].size(); j++)
                member_list[member_matrix[i][j]].push_back(i);
        //printm(member_list);
        for (unsigned long i = 0; i < member_list.size(); i++) {
            deque<long> liin;
            for (unsigned long j = 0; j < member_list[i].size(); j++) {
                compute_internal_degree_per_node(internal_degree_seq[i], member_list[i].size(), liin);
                liin.push_back(degree_seq[i] - internal_degree_seq[i]);
            }
            link_list.push_back(liin);
        }
        // now there is the check for the even node (it means that the internal degree of each group has to be even and we want to assure that, otherwise the degree_seq has to change) ----------------------------
        // ------------------------ this is done to check if the sum of the internal degree is an even number. if not, the program will change it in such a way to assure that.
        for (unsigned long i = 0; i < member_matrix.size(); i++) {
            long internal_cluster = 0;
            for (unsigned long j = 0; j < member_matrix[i].size(); j++) {
                long right_index = lower_bound(member_list[member_matrix[i][j]].begin(), member_list[member_matrix[i][j]].end(), i) - member_list[member_matrix[i][j]].begin();
                internal_cluster += link_list[member_matrix[i][j]][right_index];
            }
            if (internal_cluster % 2 != 0) {
                bool default_flag = false;
                if (excess)
                    default_flag = true;
                else if (defect)
                    default_flag = false;
                else if (ran4() > 0.5)
                    default_flag = true;
                if (default_flag) {
                    // if this does not work in a reasonable time the degree sequence will be changed
                    for (unsigned long j = 0; j < member_matrix[i].size(); j++) {
                        unsigned long random_mate = member_matrix[i][irand(member_matrix[i].size() - 1)];
                        long right_index = lower_bound(member_list[random_mate].begin(), member_list[random_mate].end(), i) - member_list[random_mate].begin();
                        if ((static_cast<unsigned long>(link_list[random_mate][right_index]) < member_matrix[i].size() - 1) && (link_list[random_mate][link_list[random_mate].size() - 1] > 0)) {
                            link_list[random_mate][right_index]++;
                            link_list[random_mate][link_list[random_mate].size() - 1]--;
                            break;
                        }
                    }
                }
                else {
                    for (unsigned long j = 0; j < member_matrix[i].size(); j++) {
                        long random_mate = member_matrix[i][irand(member_matrix[i].size() - 1)];
                        long right_index = lower_bound(member_list[random_mate].begin(), member_list[random_mate].end(), i) - member_list[random_mate].begin();
                        if (link_list[random_mate][right_index] > 0) {
                            link_list[random_mate][right_index]--;
                            link_list[random_mate][link_list[random_mate].size() - 1]++;
                            break;
                        }
                    }
                }
            }
        }
        // ------------------------ this is done to check if the sum of the internal degree is an even number. if not, the program will change it in such a way to assure that.
        {
            set<long> first;
            for (long i = 0; i < num_nodes; i++)
                E.push_back(first);
        }
        for (unsigned long i = 0; i < member_matrix.size(); i++) {
            deque<long> internal_degree_i;
            for (unsigned long j = 0; j < member_matrix[i].size(); j++) {
                long right_index = lower_bound(member_list[member_matrix[i][j]].begin(), member_list[member_matrix[i][j]].end(), i) - member_list[member_matrix[i][j]].begin();
                internal_degree_i.push_back(link_list[member_matrix[i][j]][right_index]);
            }
            if (build_subgraph(E, member_matrix[i], internal_degree_i) == -1)
                return -1;
        }
        return 0;
    }
    bool they_are_mate(long a, long b, const deque<deque<long> > & member_list) {
        for (unsigned long i = 0; i < member_list[a].size(); i++) {
            if (binary_search(member_list[b].begin(), member_list[b].end(), member_list[a][i]))
                return true;
        }
        return false;
    }
    long connect_all_the_parts(deque<set<long> > & E, const deque<deque<long> > & member_list, const deque<deque<long> > & link_list) {
        deque<long> degrees;
        for (unsigned long i = 0; i < link_list.size(); i++)
            degrees.push_back(link_list[i][link_list[i].size() - 1]);
        deque<set<long> > en; // this is the en of the subgraph
        {
            set<long> first;
            for (unsigned long i = 0; i < member_list.size(); i++)
                en.push_back(first);
        }
        multimap <long, long> degree_node;
        for (unsigned long i = 0; i < degrees.size(); i++)
            degree_node.insert(degree_node.end(), make_pair(degrees[i], i));
        long var = 0;
        while (degree_node.size() > 0) {
            multimap<long, long>::iterator itlast = degree_node.end();
            itlast--;
            multimap <long, long>::iterator itit = itlast;
            deque <multimap<long, long>::iterator> erasenda;
            long inserted = 0;
            for (long i = 0; i < itlast->first; i++) {
                if (itit != degree_node.begin()) {
                    itit--;
                    en[itlast->second].insert(itit->second);
                    en[itit->second].insert(itlast->second);
                    inserted++;
                    erasenda.push_back(itit);
                }
                else
                    break;
            }
            for (unsigned long i = 0; i < erasenda.size(); i++) {
                if (erasenda[i]->first > 1)
                    degree_node.insert(make_pair(erasenda[i]->first - 1, erasenda[i]->second));
                degree_node.erase(erasenda[i]);
            }
            var += itlast->first - inserted;
            degree_node.erase(itlast);
        }
        // this is to randomize the subgraph -------------------------------------------------------------------
        for (unsigned long node_a = 0; node_a < degrees.size(); node_a++) for (unsigned long krm = 0; krm < en[node_a].size(); krm++) {
                unsigned long random_mate = irand(degrees.size() - 1);
                while (random_mate == node_a)
                    random_mate = irand(degrees.size() - 1);
                if (en[node_a].insert(random_mate).second) {
                    deque <long> out_nodes;
                    for (set<long>::iterator it_est = en[node_a].begin(); it_est != en[node_a].end(); it_est++) if ((*it_est) != random_mate)
                            out_nodes.push_back(*it_est);
                    long old_node = out_nodes[irand(out_nodes.size() - 1)];
                    en[node_a].erase(old_node);
                    en[random_mate].insert(node_a);
                    en[old_node].erase(node_a);
                    deque <long> not_common;
                    for (set<long>::iterator it_est = en[random_mate].begin(); it_est != en[random_mate].end(); it_est++)
                        if ((old_node != (*it_est)) && (en[old_node].find(*it_est) == en[old_node].end()))
                            not_common.push_back(*it_est);
                    long node_h = not_common[irand(not_common.size() - 1)];
                    en[random_mate].erase(node_h);
                    en[node_h].erase(random_mate);
                    en[node_h].insert(old_node);
                    en[old_node].insert(node_h);
                }
            }
        // now there is a rewiring process to avoid "mate nodes" (nodes with al least one membership in common) to link each other
        long var_mate = 0;
        for (unsigned long i = 0; i < degrees.size(); i++) for (set<long>::iterator itss = en[i].begin(); itss != en[i].end(); itss++) if (they_are_mate(i, *itss, member_list)) {
                    var_mate++;
                }
        //cout<<"var mate = "<<var_mate<<endl;
        long stopper_mate = 0;
        long mate_trooper = 10;
        while (var_mate > 0) {
            //cout<<"var mate = "<<var_mate<<endl;
            long best_var_mate = var_mate;
            // ************************************************  rewiring
            for (unsigned long a = 0; a < degrees.size(); a++) for (set<long>::iterator its = en[a].begin(); its != en[a].end(); its++) if (they_are_mate(a, *its, member_list)) {
                        unsigned long b = *its;
                        unsigned long stopper_m = 0;
                        while (true) {
                            stopper_m++;
                            unsigned long random_mate = irand(degrees.size() - 1);
                            while (random_mate == a || random_mate == b)
                                random_mate = irand(degrees.size() - 1);
                            if (!(they_are_mate(a, random_mate, member_list)) && (en[a].find(random_mate) == en[a].end())) {
                                deque <long> not_common;
                                for (set<long>::iterator it_est = en[random_mate].begin(); it_est != en[random_mate].end(); it_est++)
                                    if ((b != (*it_est)) && (en[b].find(*it_est) == en[b].end()))
                                        not_common.push_back(*it_est);
                                if (not_common.size() > 0) {
                                    long node_h = not_common[irand(not_common.size() - 1)];
                                    en[random_mate].erase(node_h);
                                    en[random_mate].insert(a);
                                    en[node_h].erase(random_mate);
                                    en[node_h].insert(b);
                                    en[b].erase(a);
                                    en[b].insert(node_h);
                                    en[a].insert(random_mate);
                                    en[a].erase(b);
                                    if (!they_are_mate(b, node_h, member_list))
                                        var_mate -= 2;
                                    if (they_are_mate(random_mate, node_h, member_list))
                                        var_mate -= 2;
                                    break;
                                }
                            }
                            if (stopper_m == en[a].size())
                                break;
                        }
                        break; // this break is done because if you erased some link you have to stop this loop (en[i] changed)
                    }
            // ************************************************  rewiring
            if (var_mate == best_var_mate) {
                stopper_mate++;
                if (stopper_mate == mate_trooper)
                    break;
            } else
                stopper_mate = 0;
            //cout<<"var mate = "<<var_mate<<endl;
        }
        //cout<<"var mate = "<<var_mate<<endl;
        for (unsigned long i = 0; i < en.size(); i++) {
            for (set<long>::iterator its = en[i].begin(); its != en[i].end(); its++) if (i<static_cast<unsigned long>(*its)) {
                    E[i].insert(*its);
                    E[*its].insert(i);
                }
        }
        return 0;
    }
    long internal_kin(deque<set<long> > & E, const deque<deque<long> > & member_list, long i) {
        long var_mate2 = 0;
        for (set<long>::iterator itss = E[i].begin(); itss != E[i].end(); itss++) if (they_are_mate(i, *itss, member_list))
                var_mate2++;
        return var_mate2;
    }
    long internal_kin_only_one(set<long> & E, const deque<long> & member_matrix_j) { // return the overlap between E and member_matrix_j
        long var_mate2 = 0;
        for (set<long>::iterator itss = E.begin(); itss != E.end(); itss++) {
            if (binary_search(member_matrix_j.begin(), member_matrix_j.end(), *itss))
                var_mate2++;
        }
        return var_mate2;
    }
    long erase_links(deque<set<long> > & E, const deque<deque<long> > & member_list, const bool excess, const bool defect, const double mixing_parameter) {
        long num_nodes = member_list.size();
        long eras_add_times = 0;
        if (excess) {
            for (long i = 0; i < num_nodes; i++) {
                while ((E[i].size() > 1) && double(internal_kin(E, member_list, i)) / E[i].size() < 1 - mixing_parameter) {
                    //---------------------------------------------------------------------------------
                    cout << "degree sequence changed to respect the option -sup ... " << ++eras_add_times << endl;
                    deque<long> deqar;
                    for (set<long>::iterator it_est = E[i].begin(); it_est != E[i].end(); it_est++)
                        if (!they_are_mate(i, *it_est, member_list))
                            deqar.push_back(*it_est);
                    if (deqar.size() == E[i].size()) { // this shouldn't happen...
                        cerr << "sorry, something went wrong: there is a node which does not respect the constraints. (option -sup)" << endl;
                        return -1;
                    }
                    long random_mate = deqar[irand(deqar.size() - 1)];
                    E[i].erase(random_mate);
                    E[random_mate].erase(i);
                }
            }
        }
        if (defect) {
            for (long i = 0; i < num_nodes; i++)
                while ((E[i].size() < E.size()) && double(internal_kin(E, member_list, i)) / E[i].size() > 1 - mixing_parameter) {
                    //---------------------------------------------------------------------------------
                    cout << "degree sequence changed to respect the option -inf ... " << ++eras_add_times << endl;
                    long stopper_here = num_nodes;
                    long stopper_ = 0;
                    long random_mate = irand(num_nodes - 1);
                    while (((they_are_mate(i, random_mate, member_list)) || E[i].find(random_mate) != E[i].end()) && (stopper_ < stopper_here)) {
                        random_mate = irand(num_nodes - 1);
                        stopper_++;
                    }
                    if (stopper_ == stopper_here) { // this shouldn't happen...
                        cerr << "sorry, something went wrong: there is a node which does not respect the constraints. (option -inf)" << endl;
                        return -1;
                    }
                    E[i].insert(random_mate);
                    E[random_mate].insert(i);
                }
        }
        //------------------------------------ Erasing links   ------------------------------------------------------
        return 0;
    }
    long print_network(deque<set<long> > & E, const deque<deque<long> > & member_list, const deque<deque<long> > & member_matrix, deque<long> & num_seq) {
        long edges = 0;
        long num_nodes = member_list.size();
        deque<double> double_mixing;
        for (unsigned long i = 0; i < E.size(); i++) {
            double one_minus_mu = double(internal_kin(E, member_list, i)) / E[i].size();
            double_mixing.push_back(1. - one_minus_mu);
            edges += E[i].size();
        }
        //cout<<"\n----------------------------------------------------------"<<endl;
        //cout<<endl;
        double density = 0;
        double sparsity = 0;
        for (unsigned long i = 0; i < member_matrix.size(); i++) {
            double media_long = 0;
            double media_est = 0;
            for (unsigned long j = 0; j < member_matrix[i].size(); j++) {
                double kinj = double(internal_kin_only_one(E[member_matrix[i][j]], member_matrix[i]));
                media_long += kinj;
                media_est += E[member_matrix[i][j]].size() - double(internal_kin_only_one(E[member_matrix[i][j]], member_matrix[i]));
            }
            double pair_num = (member_matrix[i].size()*(member_matrix[i].size() - 1));
            double pair_num_e = ((num_nodes - member_matrix[i].size())*(member_matrix[i].size()));
            if (pair_num != 0)
                density += media_long / pair_num;
            if (pair_num_e != 0)
                sparsity += media_est / pair_num_e;
        }
        density = density / member_matrix.size();
        sparsity = sparsity / member_matrix.size();
        ofstream out1("network.dat");
        for (unsigned long u = 0; u < E.size(); u++) {
            set<long>::iterator itb = E[u].begin();
            while (itb != E[u].end())
                out1 << u + 1 << "\t" << *(itb++) + 1 << endl;
        }
        ofstream out2("community.dat");
        for (unsigned long i = 0; i < member_list.size(); i++) {
            out2 << i + 1 << "\t";
            for (unsigned long j = 0; j < member_list[i].size(); j++)
                out2 << member_list[i][j] + 1 << " ";
            out2 << endl;
        }
        cout << "\n\n---------------------------------------------------------------------------" << endl;
        cout << "network of " << num_nodes << " vertices and " << edges / 2 << " edges" << ";\t average degree = " << double(edges) / num_nodes << endl;
        cout << "\naverage mixing parameter: " << average_func(double_mixing) << " +/- " << sqrt(variance_func(double_mixing)) << endl;
        cout << "p_in: " << density << "\tp_out: " << sparsity << endl;
        ofstream statout("statistics.dat");
        deque<long> degree_seq;
        for (unsigned long i = 0; i < E.size(); i++)
            degree_seq.push_back(E[i].size());
        statout << "degree distribution (probability density function of the degree in logarithmic bins) " << endl;
        log_histogram(degree_seq, statout, 10);
        statout << "\ndegree distribution (degree-occurrences) " << endl;
        int_histogram(degree_seq, statout);
        statout << endl << "--------------------------------------" << endl;
        statout << "community distribution (size-occurrences)" << endl;
        int_histogram(num_seq, statout);
        statout << endl << "--------------------------------------" << endl;
        statout << "mixing parameter" << endl;
        not_norm_histogram(double_mixing, statout, 20, 0, 0);
        statout << endl << "--------------------------------------" << endl;
        cout << endl << endl;
        return 0;
    }
    long benchmark(bool excess, bool defect, long num_nodes, double average_k, long max_degree, double tau, double tau2, double mixing_parameter, long overlapping_nodes, long overlap_membership, long nmin, long nmax, bool fixed_range, deque<set<long> >& E, deque<deque<long> >& member_list, deque<deque<long> >& link_list) {
        double dmin = solve_dmin(max_degree, average_k, -tau);
        if (dmin == -1)return -1;
        long min_degree = int(dmin);
        double media1 = integer_average(max_degree, min_degree, tau);
        double media2 = integer_average(max_degree, min_degree + 1, tau);
        if (fabs(media1 - average_k) > fabs(media2 - average_k)) min_degree++;
        // range for the community sizes
        if (!fixed_range) {
            nmax = max_degree;
            nmin = max(int(min_degree), 3);
        }
        deque <long> degree_seq; //  degree sequence of the nodes
        deque <double> cumulative;
        powerlaw(max_degree, min_degree, tau, cumulative);
        for (long i = 0; i < num_nodes; i++) {
            long nn = lower_bound(cumulative.begin(), cumulative.end(), ran4()) - cumulative.begin() + min_degree;
            degree_seq.push_back(nn);
        }
        sort(degree_seq.begin(), degree_seq.end());
        if (deque_int_sum(degree_seq) % 2 != 0) degree_seq[max_element(degree_seq.begin(), degree_seq.end()) - degree_seq.begin()]--;
        deque<deque<long> > member_matrix;
        deque<long> num_seq;
        deque<long> internal_degree_seq;
        if (internal_degree_and_membership(mixing_parameter, overlapping_nodes, overlap_membership, num_nodes, member_matrix, excess, defect, degree_seq, num_seq, internal_degree_seq, fixed_range, nmin, nmax, tau2) == -1) return -1;
        //deque<set<long> > E;					// E is the adjacency matrix written in form of list of edges
        //deque<deque<long> > member_list;		// row i cointains the memberships of node i
        //deque<deque<long> > link_list;		// row i cointains degree of the node i respect to member_list[i][j]; there is one more number that is the external degree
        if (build_subgraphs(E, member_matrix, member_list, link_list, internal_degree_seq, degree_seq, excess, defect) == -1) return -1;
        connect_all_the_parts(E, member_list, link_list);
        //print_network(E, member_list, member_matrix, num_seq);
        return 0;
    }
};
#endif