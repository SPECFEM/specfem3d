//
//    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
//                    Simon Stähler, Kasra Hosseini, Stefanie Hempel
//
//    This file is part of AxiSEM.
//    It is distributed from the webpage <http://www.axisem.info>
//
//    AxiSEM is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    AxiSEM is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with AxiSEM.  If not, see <http://www.gnu.org/licenses/>.
//

#include <stdio.h>
#include <pthread.h>
#include <unistd.h>

// compare https://computing.llnl.gov/tutorials/pthreads/

// careful with global variable, this only works if the calling programm is not threaded!
pthread_t thread;

// stub - so c knows what the function looks like
extern void __nc_routines_MOD_nc_dump_strain_to_disk();

// tread that calls the IO function from fortran
void *cfunc_thread(void* valp)
{
   int val;
   val = *((int*)valp);
   nc_dump_strain_to_disk();
   pthread_exit(NULL);
}

// create IO thread, to be called from fortran
void c_spawn_dumpthread(int* val){
   pthread_attr_t attr; 
   pthread_attr_init(&attr);
   pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
   pthread_create(&thread, &attr, cfunc_thread, (void *)val);
}

// wait for the IO thread to finish, to be called from fortran
// global thead variable allows to come back to the thread
void c_wait_for_io() {
   void *status;
   pthread_join(thread, &status);
}

