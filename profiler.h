/* 
 * File:   profiler.h
 * Author: bharath
 *
 * Created on May 11, 2012, 10:11 AM
 */

#ifndef PROFILER_H
#define	PROFILER_H

#include "typedefs.h"

namespace CDLib
{
        class timer_rt
    {
    private:
        timespec start;
        timespec stop;
        double all_time;
    public:
        timer_rt(): start(),stop(),all_time(0.0){}
        void start_clock() { clock_gettime(CLOCK_MONOTONIC, &start); }
        double run_time() const { return (stop.tv_sec - start.tv_sec) + (double)(stop.tv_nsec - start.tv_nsec) / (double)BILLION; }
        void stop_clock() 
        { 
            clock_gettime(CLOCK_MONOTONIC, &stop); 
            all_time += run_time();
        }
        double total_time() const { return all_time; }
    };
    
    id_type get_mem_usage()
    {
        rusage usage; 
        getrusage(RUSAGE_SELF,&usage);
        return usage.ru_maxrss;
    }
};



#endif	/* PROFILER_H */

