/* Watchdog routine written by David Radice for the Einstein Toolkit
   It checks whether a process has made progress within a predefined time.
   If not, it will abort the program.
   License LGPL.
*/
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include <pthread.h>

// If a process does not progress after this time, the program will abort
static const int checkevery=3600;

static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
static time_t timestamp;
static struct {
    int timeout_sec;
    int mpi_rank;
} param;

static pthread_t dog;
static void * patrol(void * arg) {
    time_t time_old, time_new, ltime;
    char tstamp[128];

    pthread_mutex_lock(&mutex);
    time_old = timestamp;
    pthread_mutex_unlock(&mutex);

    if(0 == param.mpi_rank) {
        ltime = time(NULL);
        asctime_r(localtime(&ltime), &tstamp[0]);
        tstamp[24] = '\0';
        fprintf(stderr, "[WATCHDOG (%s)] Starting.\n", tstamp);
        fflush(stderr);
    }

    while(true) {
        unsigned int left = param.timeout_sec;
        while(left > 0) {
            left = sleep(left);
        }

        pthread_mutex_lock(&mutex);
        time_new = timestamp;
        pthread_mutex_unlock(&mutex);

        ltime = time(NULL);
        asctime_r(localtime(&ltime), &tstamp[0]);
        tstamp[24] = '\0';
        if(time_new == time_old) {
            fprintf(stderr, "[WATCHDOG (%s)] Rank %d is not progressing.\n",
                    tstamp, param.mpi_rank);
            fprintf(stderr, "[WATCHDOG (%s)] Terminating...\n", tstamp);
            fflush(stderr);
            abort();
        }
        else {
            if(0 == param.mpi_rank) {
                fprintf(stderr, "[WATCHDOG (%s)] Everything is fine.\n",
                        tstamp);
                fflush(stderr);
            }
            time_old = time_new;
        }
    }
}

void WatchDog_update() {

    pthread_mutex_lock(&mutex);
    timestamp = time(NULL);
    pthread_mutex_unlock(&mutex);
}

void WatchDog_create(int myrank){ //, int checkevery) {

    WatchDog_update();

        param.timeout_sec = checkevery;
        param.mpi_rank = myrank;
        int ierr = pthread_create(&dog, NULL, patrol, NULL);
        if(ierr) {
            fprintf(stderr, "Dispatching watch dog thread failed");
            fflush(stderr);
            abort();
        }
}
