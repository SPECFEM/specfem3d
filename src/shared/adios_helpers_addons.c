/*
!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
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
!=====================================================================
*/

/* ----------------------------------------------------------------------------- */

// ADIOS2 add-on wrappers

/* ----------------------------------------------------------------------------- */
#if defined(USE_ADIOS2)

// adios2 provides some of the functionality only to C bindings, but not fortran.
// here we add some wrappers for such function calls to be able to call them from our fortran routines.

#include <stdio.h>
#include <stdlib.h>
#include "config.h"

#include "adios2_c.h"

void
FC_FUNC_(get_adios2_all_variables_count,GET_ADIOS2_ALL_VARIABLES_COUNT)(int *vars_count, adios2_io **adios_group) {

  adios2_variable **adios_vars;
  adios2_error ret;
  size_t nvars;

  ret = adios2_inquire_all_variables(&adios_vars, &nvars, *adios_group);
  if (ret != 0){ fprintf(stderr, "Error adios2 inquire all variables\n"); exit(-1);}

  *vars_count = (int) nvars;
  free(adios_vars);
}

/* ------------------------------------------------------------- */

void
FC_FUNC_(show_adios2_all_variables,SHOW_ADIOS2_ALL_VARIABLES)(adios2_io **adios_group) {

  adios2_variable **adios_vars;
  adios2_error ret;
  size_t nvars;

  ret = adios2_inquire_all_variables(&adios_vars, &nvars, *adios_group);
  if (ret != 0){ fprintf(stderr, "Error adios2 inquire all variables\n"); exit(-1);}

  // user output
  printf("variables: %i\n",(int)nvars);
  for (int i = 0; i < nvars; i++){
    adios2_variable *var = adios_vars[i];
    if (var == NULL) continue;

    size_t len;
    ret = adios2_variable_name(NULL, &len, var);
    if (ret != 0){ fprintf(stderr, "Error adios2 getting variable name length\n"); exit(-1);}

    if (len > 0){
      char name[len+1];
      ret = adios2_variable_name(name, &len, var);
      if (ret != 0){ fprintf(stderr, "Error adios2 getting variable name\n"); exit(-1);}
      name[len] = '\0';
      printf("  id %i : name = %s\n",i,name);
    }else{
      printf("  id %i : has zero name length\n",i);
    }
  }
}

/* ------------------------------------------------------------- */

void
FC_FUNC_(get_adios2_all_attributes_count,GET_ADIOS2_ALL_ATTRIBUTES_COUNT)(int *attrs_count, adios2_io **adios_group) {

  adios2_attribute **adios_attrs;
  adios2_error ret;
  size_t nattrs;

  ret = adios2_inquire_all_attributes(&adios_attrs, &nattrs, *adios_group);
  if (ret != 0){ fprintf(stderr, "Error adios2 inquire all variables\n"); exit(-1);}

  *attrs_count = (int) nattrs;
  free(adios_attrs);
}

/* ------------------------------------------------------------- */

void
FC_FUNC_(show_adios2_all_attributes,SHOW_ADIOS2_ALL_ATTRIBUTES)(adios2_io **adios_group) {

  adios2_attribute **adios_attrs;
  adios2_error ret;
  size_t nattrs;

  ret = adios2_inquire_all_attributes(&adios_attrs, &nattrs, *adios_group);
  if (ret != 0){ fprintf(stderr, "Error adios2 inquire all variables\n"); exit(-1);}

  // user output
  printf("attributes: %i\n",(int)nattrs);
  for (int i = 0; i < nattrs; i++){
    adios2_attribute *attr = adios_attrs[i];

    size_t len;
    ret = adios2_attribute_name(NULL, &len, attr);
    if (ret != 0){ fprintf(stderr, "Error adios2 getting attribute name length\n"); exit(-1);}

    if (len > 0){
      char name[len+1];
      ret = adios2_attribute_name(name, &len, attr);
      if (ret != 0){ fprintf(stderr, "Error adios2 getting attribute name\n"); exit(-1);}
      name[len] = '\0';
      printf("  id %i : name = %s\n",i,name);
    }else{
      printf("  id %i : has zero name length\n",i);
    }
  }

}

#endif

