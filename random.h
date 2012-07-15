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
    template <typename T>
    class RandomGenerator
    {
        private:
            T t_max;
            T t_min;

        public:
            static void seed(unsigned long seed) { srand(seed); }
            RandomGenerator() :t_max(1),t_min(0) { }
            RandomGenerator(T min, T max) :t_max(max),t_min(min) { }
            RandomGenerator(T min, T max, bool autoseed) :t_max(max),t_min(min) { if(autoseed) seed(time(NULL)); }
            T next()  { return (T)((t_max - t_min)*((double)rand()/RAND_MAX) + t_min); }
            T exp_next(double alpha, double beta) { return (T)(beta*pow(next(),-alpha)); }
    };

    class RandomStringGenerator
    {
    private:
        RandomGenerator<size_t> rgl_size;
        RandomGenerator<char> rgc_char;
    public:
        RandomStringGenerator() : rgl_size(1,2,1), rgc_char('a','z'+1,1) { }
        RandomStringGenerator(size_t min_size,size_t max_size,char start_char, char end_char)  : rgl_size(min_size,max_size,1), rgc_char(start_char,end_char+1,1) { }
        void next(string& outStr)
        {
            size_t length = rgl_size.next();
            for(size_t i = 0;i<length;i++)
                outStr.push_back(rgc_char.next());
        }
    };
//<<<<<<< .mine
    
    
    
    /*UniformRandomGenerator is used to generate random number uniformly distributed over certain range.
     if you want to generate # between 0 & 1 then use double as a typename.
     if you want to generate # between 0 & 'any natural #' then use any integer format as a typename.
     this function can generate different sequence every time. 
     declare UniformRandomGenerator class object out side all loops. */
//    
    template <typename T1>
    class UniformRandomGenerator
    {
    private:
        T1 t_max;
        T1 t_min;
    public:
        UniformRandomGenerator() {srand((unsigned)time(0));}
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
        UniformRandomGenerator<double> randdouble;
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
