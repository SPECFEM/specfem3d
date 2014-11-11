/*
 * Copyright [2014] [Matthieu Lefebvre]
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>

#include "sep_header.h"
#include "parse_sep.h"
#include "check_errors.h"

#include "convert_topo.h"

/*****************************************************************************/

void convert_topo_sep_to_ascii(char *sep_name, char *ascii_name) {
  /* Get sep header from sep file 'sep_name' */
  sep_header_t sep_header;
  parse_sep_vs_header(sep_name, &sep_header);
  show_header(sep_header);
  /* Read topo value from sep binary file */
  int *sep_topo = (int *) malloc(sep_header.n1 * sep_header.n2
                                    * sizeof(int));
  read_topo_sep(&sep_header, sep_topo);
  /* Create ascii topo according to geometry */
  float *ascii_topo = (float *) malloc(sep_header.n1 * sep_header.n2
                                      * sizeof (float));
  /* sep topo index to ascii topo depth (m) */
  /* Interpolation is done in meshfem 3d*/
  assign_ascii_topo_values(&sep_header, sep_topo, ascii_topo);
  write_ascii_topo(ascii_name, ascii_topo, sep_header.n1 * sep_header.n2);
  /* Clean */
  free(sep_topo);
  free(ascii_topo);
}

void read_topo_sep(sep_header_t *header, int *sep_topo) {
  FILE *fd = fopen(header->in, "r");
  if (!fd) EXIT_ON_ERR("fopen");

  fread(sep_topo, sizeof(int), header->n1 * header->n2, fd);

  fclose(fd);
}

void assign_ascii_topo_values(sep_header_t *header,
                              int *sep_topo, float *ascii_topo) {
  int nx = header->n1;
  int ny = header->n2;

  float oz = header->o3;
  float dz = header->d3;

  for (int i = 0; i < nx*ny; ++i) {
    ascii_topo[i] = - dz*(oz + sep_topo[i]);
  }
}

void write_ascii_topo(char *ascii_name, float *ascii_topo, int n) {
  FILE *fd = fopen(ascii_name, "w");
  if (fd == NULL)
    EXIT_ON_ERR("fopen");

  for (int i = 0; i < n; ++i) {
    fprintf(fd, "%f\n", ascii_topo[i]);
  }

  fclose(fd);
}
