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
#include <getopt.h>
#include <string.h>

#include "sep_header.h"
#include "parse_sep.h"
#include "convert_topo.h"

/*****************************************************************************/
/**
  * \brief Parse arguments form the command line.
  *
  * \param argc
  * \param argv
  * \param vs_name From '-i' option.
  * \param topo_name From '-o' option.
  *
  * \return Success if all args are correct.
  */
int parse_args(int argc, char *argv[], char **vs_name, char **topo_name);

/*****************************************************************************/
int main(int argc, char *argv[]) {
  char *vs_name, *topo_name;
  parse_args(argc, argv, &vs_name, &topo_name);

  sep_header_t vs_header;
  parse_sep_vs_header(vs_name, &vs_header);
  show_header(vs_header);

  int *topo = (int*) malloc(vs_header.n1 * vs_header.n2 * sizeof(int));
  parse_sep_vs_binary(&vs_header, topo);

  /* Create ascii topo according to geometry */
  float *ascii_topo = (float *) malloc(vs_header.n1 * vs_header.n2
                                      * sizeof (float));
  /* sep topo index to ascii topo depth (m) */
  /* Interpolation is done in meshfem 3d*/
  assign_ascii_topo_values(&vs_header, topo, ascii_topo);
  write_ascii_topo(topo_name, ascii_topo, vs_header.n1 * vs_header.n2);
  /* Clean */
  free(topo);
  free(ascii_topo);

  if (vs_name) free(vs_name);
  if (topo_name) free(topo_name);

  return EXIT_SUCCESS;
}

/*****************************************************************************/
int parse_args(int argc, char *argv[], char **vs_name, char **topo_name) {
  int choice;
  while (1) {
    static struct option long_options[] =
    {
      /* Argument styles: no_argument, required_argument, optional_argument */
      {"help", optional_argument, 0, 'h'},
      {"input", required_argument, 0, 'i'},
      {"output",required_argument, 0, 'o'},
      {0,0,0,0}
    };

    int option_index = 0;

    /* Argument parameters:
      no_argument: " "
      required_argument: ":"
      optional_argument: "::" */

    choice = getopt_long( argc, argv, "h::i:o:",
          long_options, &option_index);

    if (choice == -1)
      break;

    switch( choice ) {
      case 'v':
        break;
      case 'i':
        *vs_name = malloc (strlen(optarg) * sizeof(char));
        strncpy(*vs_name, optarg, strlen(optarg));
        break;
      case 'o':
        *topo_name = malloc (strlen(optarg) * sizeof(char));
        strncpy(*topo_name, optarg, strlen(optarg));
        break;
      /* For all following cases, print a hel string. */
      case 'h':
      case '?':
      default:
        printf("Usage: ./bin/create_ascii_bathy_from_sep_vs "
               "-i [input_sep_header.H] -o [output_bathy.txt]\n");
        exit(EXIT_FAILURE);
    }
  }
  /* Deal with non-option arguments here */
  if ( optind < argc ) {
    while ( optind < argc ) {
      continue;
    }
  }
  return EXIT_SUCCESS;
}
