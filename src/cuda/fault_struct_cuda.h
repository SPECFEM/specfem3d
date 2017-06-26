/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  3 . 0
 !               ---------------------------------------
 !
 !     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
 !                        Princeton University, USA
 !                and CNRS / University of Marseille, France
 !                 (there are currently many more authors!)
 ! (c) Princeton University and CNRS / University of Marseille, July 2012
 !
 ! This program is free software; you can redistribute it and/or modify
 ! it under the terms of the GNU General Public License as published by
 ! the Free Software Foundation; either version 3 of the License, or
 ! (at your option) any later version.
 !
 ! This program is distributed in the hope that it will be useful,
 ! but WITHOUT ANY WARRANTY; without even the implied warranty of
 ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ! GNU General Public License for more details.
 !
 ! You should have received a copy of the GNU General Public License along
 ! with this program; if not, write to the Free Software Foundation, Inc.,
 ! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 !
 !=====================================================================*/
/* ----------------------------------------------------------------------------------------------- */

// mesh pointer wrapper structure

/* ----------------------------------------------------------------------------------------------- */

#include "stdlib.h"
#include "config.h"

typedef float realw;

typedef struct swf_type_{
  int kind;
  int healing ;
  realw *Dc,*mus,*mud,*theta, *T, *Coh;
}Swf_type;

typedef struct rsf_type_{
  int StateLaw ; // By default using aging law
  realw *V0,*f0,*L,*V_init,*a,*b,*theta,*T,*C,*fw,*Vw;
}Rsf_type;

typedef struct fault_data_{
  int NRECORD;
  int NT;
  int* iglob; //record the global index of the station
  realw* dataT; //data pointer to the array
}Fault_data;

typedef struct fault_ {

  // mesh resolution
  int NSPEC_AB;
  int NGLOB_AB;

  Fault_data* output_dataT;
  realw *T0, *T, *B, *V, *D;
  realw *R, *invM1, *invM2, *Z ;
  int *ibulk1,*ibulk2;
  Swf_type swf;
  Rsf_type rsf;

}Fault ;

typedef struct fault_solver_dynamics_{
  Fault* faults;
  realw v_healing;
  realw v_rupt;
  int NTOUT,NSNAP;
  int Nbfaults;
  int RATE_AND_STATE;
}Fault_solver_dynamics;

enum Simulationtype
{
  Slipweakening=0,
  RateandstateAging=1,
  RateandstateStron=2
};
