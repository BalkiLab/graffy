/*
 * File:   profiler.h
 * Author: bharath
 *
 * Created on May 11, 2012, 10:11 AM
 */

#ifndef PROFILER_H
#define	PROFILER_H

#include "typedefs.h"

namespace CDLib {

    class timer_rt {
    private:
        timespec start;
        timespec stop;
        time_t clock_start, clock_stop;
        double runtime, all_time;

        string get_time(const time_t * st_time) {
            struct tm *read_time = localtime(st_time);
            ostringstream oss;
            oss << read_time->tm_mday << "/" << (1 + read_time->tm_mon) << "/" << (1900 + read_time->tm_year) << " " << read_time->tm_hour << ":" << read_time->tm_min << ":" << read_time->tm_sec;
            return oss.str();
        }

        inline double diff(timespec& end, timespec& begin) {
            return (end.tv_sec - begin.tv_sec) + (double) (end.tv_nsec - begin.tv_nsec) / (double) BILLION;
        }
    public:

        timer_rt() : start(), stop(), clock_start(), clock_stop(), all_time(0.0) {
        }

        void start_clock() {
            clock_start = time(NULL);
            clock_gettime(CLOCK_MONOTONIC, &start);
        }

        double run_time() const {
            return runtime;
        }

        double time_elapsed() {
            timespec cstop;
            clock_gettime(CLOCK_MONOTONIC, &cstop);
            return diff(cstop, start);
        }

        void stop_clock() {
            clock_stop = time(NULL);
            clock_gettime(CLOCK_MONOTONIC, &stop);
            runtime = diff(stop, start);
            all_time += runtime;
        }

        double total_time() const {
            return all_time;
        }

        string start_time() {
            return get_time(&clock_start);
        }

        string stop_time() {
            return get_time(&clock_stop);
        }
    };

    id_type get_mem_usage() {
        rusage usage;
        getrusage(RUSAGE_SELF, &usage);
        return usage.ru_maxrss;
    }
};



#endif	/* PROFILER_H */

