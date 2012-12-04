/* 
 * File:   random.h
 * Author: bharath
 *
 * Created on 16 April, 2012, 3:53 PM
 */

#ifndef RANDOM_H
#define	RANDOM_H

#include "typedefs.h"

using namespace std;

namespace CDLib {
    
    class Uniform01RandomGeneratorDefault{
    private:
        double convert_to_01(long val) { return static_cast<double>(val)/static_cast<double>(RAND_MAX);}
    public:
        Uniform01RandomGeneratorDefault() { srand(time(NULL));}
        Uniform01RandomGeneratorDefault(long s) { srand(s); }
        inline double next() { return convert_to_01(rand()); }
    };
    
    class Uniform01RandomGeneratorMT
    {
    private:
        mt19937 gen;
    public:
        Uniform01RandomGeneratorMT() : gen(time(NULL)) {}
        Uniform01RandomGeneratorMT(long s) : gen(s) {}
        inline double next() { return static_cast<double>(static_cast<double>(gen())/static_cast<double>((unsigned int)4294967295)); }
    };
    
    template <typename T,class U01RG>
    class UniformRandomGenerator
    {
    private:
        T t_max;
        T t_min;
        U01RG gen;
    public:
            UniformRandomGenerator(long seed_value) :t_max(1),t_min(0),gen(seed_value){ }
            UniformRandomGenerator() :t_max(1),t_min(0),gen(srand(time(NULL))){ }
            UniformRandomGenerator(long seed_value,T min, T max) :t_max(max),t_min(min),gen(seed_value){ }
            UniformRandomGenerator(T min, T max) :t_max(max),t_min(min),gen(srand(time(NULL))){ }
            inline T next() { return static_cast<T>((t_max - t_min)*(gen.next()) + t_min);  }
    };
    
    template<class RandomAccessIterator, class U01RG>
    void UniformShuffler(RandomAccessIterator begin,RandomAccessIterator end,long seed) {
        U01RG gen(seed);
	for (id_type siz=end-begin; siz>=0; siz--) 
            swap(*(begin + static_cast<id_type>(siz*gen.next())),*(begin + siz-1));
    }
    
    template<class RandomAccessIterator>
    void UniformShufflerDefault(RandomAccessIterator begin,RandomAccessIterator end,long seed) {
        srand(seed);
        random_shuffle(begin,end);
    }
    
    template<class RandomAccessIterator>
    void UniformShufflerMT(RandomAccessIterator begin,RandomAccessIterator end,long seed) {
        mt19937 gen(seed);
        shuffle(begin,end,gen);
    }
    
    const long unlikely = 214741;
    const double R2_EPS = 1.2e-7;
    const double R2_RNMX = (1.0-R2_EPS);
    const long R2_IM1 = 2147483563;
    const long R2_IM2 = 2147483399;
    const long R2_IA1 = 40014;
    const long R2_IA2 = 40692;
    const long R2_IQ1 = 53668;
    const long R2_IQ2 = 52774;
    const long R2_IR1 = 12211;
    const long R2_IR2 = 3791;
    const long R2_NTAB = 32;
    const double R2_AM = (1.0/R2_IM1);
    const long R2_IMM1 = (R2_IM1-1);
    const long R2_NDIV = (1+R2_IMM1/R2_NTAB);
    
    class Uniform01RandomGeneratorLFR
    {
    private:
        long seed;
        long idum2;
        long iy;
        long iv[R2_NTAB];
        double r01(long *idum) {
            int j;
            long k;
            double temp;
            if(*idum<=0 || !iy){
                    if(-(*idum)<1) *idum=1*(*idum);
                    else *idum=-(*idum);
                    idum2=(*idum);
                    for(j=R2_NTAB+7;j>=0;j--){
                            k=(*idum)/R2_IQ1;
                            *idum=R2_IA1*(*idum-k*R2_IQ1)-k*R2_IR1;
                            if(*idum<0) *idum+=R2_IM1;
                            if(j<R2_NTAB) iv[j]=*idum;
                    }
                    iy=iv[0];
            }
            k=(*idum)/R2_IQ1;
            *idum=R2_IA1*(*idum-k*R2_IQ1)-k*R2_IR1;
            if(*idum<0) *idum+=R2_IM1;
            k=(idum2)/R2_IQ2;
            idum2=R2_IA2*(idum2-k*R2_IQ2)-k*R2_IR2;
            if (idum2 < 0) idum2 += R2_IM2;
            j=iy/R2_NDIV;
            iy=iv[j]-idum2;
            iv[j]=*idum;
            if(iy<1) iy+=R2_IMM1;
            if((temp=R2_AM*iy)>R2_RNMX) return R2_RNMX;
            else return temp;
        }
        
    public:
        inline void init(){
            idum2=123456789;
            iy=0;
        }
        Uniform01RandomGeneratorLFR() : seed(time(NULL)) { init();}
        Uniform01RandomGeneratorLFR(long s) : seed(s) { init();}
        inline double next() { return r01(&seed); }
    };
        
    template <class U01RG>
    class DiscretePowerLawGenerator
    {
    private:
        vector<double> cdf;
        U01RG gen;
        id_type t_min;
        id_type t_avg;
        id_type t_max;
        double exponent;
        inline double integral(double a, double b) {
            if (fabs(a+1.0)>1e-10) return (1.0/(a+1.0)*pow(b, a+1.0));
            else return (log(b));
        }
        
        inline double average(double dmax, double dmin, double gamma) {
            return (1.0/(integral(gamma, dmax)-integral(gamma, dmin)))*(integral(gamma+1, dmax)-integral(gamma+1, dmin));
        }
        
        double integer_average (id_type n, id_type min, double tau) {
            double a=0;
            for (double h=min; h<n+1; h++)
                a+= pow((1./h),tau);
            double pf=0;
            for(double i=min; i<n+1; i++) 
                pf+=1/a*pow((1./(i)),tau)*i;
            return pf;
        }
        
        void set_min(double threshold) {
            double dmax = t_max,dmed = t_avg,gamma = -exponent;
            double dmin_l=1,dmin_r=dmax;
            double average_k1=average(dmin_r, dmin_l, gamma),average_k2=dmin_r;
            if ((average_k1-dmed<0) && (average_k2-dmed>0)){
                while (fabs(average_k1-dmed)>threshold) {
                        double temp=average(dmax, ((dmin_r+dmin_l)/2.), gamma);
                        if ((temp-dmed)*(average_k2-dmed)>0) {
                                average_k2=temp;
                                dmin_r=((dmin_r+dmin_l)/2.);
                        }
                        else {
                            average_k1=temp;
                            t_min=((dmin_r+dmin_l)/2.);
                        }
                }
                t_min = static_cast<id_type>(dmin_l);
                double media1=integer_average(t_max, t_min, exponent);
                double media2=integer_average(t_max, t_min+1, exponent);
                if (fabs(media1-t_avg)>fabs(media2-t_avg))t_min++;
		
            }
        }
        
        void generate_cdf() {
            //set_min(1e-5);
            t_min = 1;
            cdf.clear();
            double a=0;		
            for (double h=t_min; h<t_max+1; h++) a+= pow((1.0/h),exponent);
            double pf=0;
            for(double i=t_min; i<t_max; i++) {
                    pf+=1/a*pow((1.0/(i)),exponent);
                    cdf.push_back(pf);
            }
        }
        
    public:
        DiscretePowerLawGenerator(id_type avg_val, id_type max_val,double exp_value,long seed_value) : cdf(),gen(seed_value),t_min(1),t_avg(avg_val),t_max(max_val),exponent(exp_value) { generate_cdf(); }
        DiscretePowerLawGenerator(id_type avg_val, id_type max_val,double exp_value) : cdf(),gen(time(NULL)),t_min(1),t_avg(avg_val),t_max(max_val),exponent(exp_value) { generate_cdf(); }
        id_type next() { return lower_bound(cdf.begin(), cdf.end(), gen.next())-cdf.begin()+t_min; }
    };
    
    template <typename T>
    class RandomGenerator
    {
        private:
            T t_max;
            T t_min;

        public:
            static void seed(unsigned long seed) 
            {
                srand(seed); 
            }
            RandomGenerator() :t_max(1),t_min(0) { }
            RandomGenerator(T min, T max) :t_max(max),t_min(min){ }
            RandomGenerator(T min, T max, bool autoseed) :t_max(max),t_min(min) { if(autoseed) seed(time(NULL)); }
            T next() { return static_cast<T>((t_max - t_min)*(static_cast<double>(rand())/RAND_MAX) + t_min);  }
            T exp_next(double alpha, double beta) { return (T)(beta*pow(next(),-alpha)); }
    };

    class RandomStringGenerator
    {
    private:
        RandomGenerator<size_t> rgl_size;
        RandomGenerator<char> rgc_char;
    public:
        RandomStringGenerator() : rgl_size(1,2,1), rgc_char('a','z'+1,(bool)1) { }
        RandomStringGenerator(size_t min_size,size_t max_size,char start_char, char end_char)  : rgl_size(min_size,max_size,1), rgc_char(start_char,end_char+1,1) { }
        void next(string& outStr)
        {
            size_t length = rgl_size.next();
            for(size_t i = 0;i<length;i++)
                outStr.push_back(rgc_char.next());
        }
    };
    
    
    /*UniformRandomGeneratorAkash is used to generate random number uniformly distributed over certain range.
     if you want to generate # between 0 & 1 then use double as a typename.
     if you want to generate # between 0 & 'any natural #' then use any integer format as a typename.
     this function can generate different sequence every time. 
     declare UniformRandomGeneratorAkash class object out side all loops. */
//    
    template <typename T1>
    class UniformRandomGeneratorAkash
    {
    private:
        T1 t_max;
        T1 t_min;
    public:
        UniformRandomGeneratorAkash() {srand((unsigned)time(0));}
        T1 next(T1 max) 
        {
            t_max=max;
            return (T1)(((double)rand()/RAND_MAX) * t_max);
        }
        T1 next(T1 min,T1 max)
        {
            t_min=min;
            t_max=max;
            return (T1)((t_max - t_min)*((double)rand()/RAND_MAX) + t_min);
        }
             
    };
    
    template <typename T1>
    class GaussianRandomGenerator
    {
    private:
        T1 Mean;
        T1 Variance;
        UniformRandomGeneratorAkash<double> randdouble;
    public:
        T1 next(T1 mean,T1 variance)
        {
            Mean=mean;
            Variance=variance;
            double u1,u2;
            u1=randdouble.next(1);
            u2=randdouble.next(1);
            double theta,R;
            
            theta = 2*3.14159265*u2;
            R = sqrt((-2)*log(u1));
                      
            double gaussianrandom_0_1;
            
            gaussianrandom_0_1 = R*cos(theta);
                        
            return (T1)(gaussianrandom_0_1*Variance) + Mean;
        }
    };
};

#endif	/* RANDOM_H */
