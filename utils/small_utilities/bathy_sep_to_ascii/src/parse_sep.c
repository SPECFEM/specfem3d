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
#include <string.h>
#include <regex.h>

#include "parse_sep.h"
#include "check_errors.h"

/*****************************************************************************/
void parse_sep_vs_header(char *vs_name, sep_header_t *sep_header) {
  char *line = NULL;
  size_t len = 0;

  FILE *fd = fopen(vs_name, "r");
  if (fd == NULL)
    EXIT_ON_ERR("fopen");

  regex_t re;
  if(regcomp(&re, "^[ \t]*\\(.*\\)[ \t]*=[ \t]*\\(.*\\)[ \t]*",
             REG_ICASE | REG_NEWLINE))
    EXIT_ON_ERR("regcomp");

  while(getline(&line, &len, fd) != -1) {

    char name[512], value[512];
    size_t nmatch = 3;
    regmatch_t pmatch[3];
    if(REG_NOMATCH == regexec(&re, line, nmatch, pmatch, 0)) {
      continue;
    } else {
      size_t name_len = pmatch[1].rm_eo - pmatch[1].rm_so;
      strncpy(name, &line[pmatch[1].rm_so], name_len);
      name[name_len] = '\0';

      size_t value_len = pmatch[2].rm_eo - pmatch[2].rm_so;
      strncpy(value, &line[pmatch[2].rm_so], pmatch[2].rm_eo - pmatch[2].rm_so);
      value[value_len] = '\0';

      if (!strcmp(name, "n1")) {
        sep_header->n1 = atoi(value);
      } else if (!strcmp(name, "n2")) {
        sep_header->n2 = atoi(value);
      } else if (!strcmp(name, "n3")) {
        sep_header->n3 = atoi(value);
      } else if (!strcmp(name, "o1")) {
        sep_header->o1 = atoi(value);
      } else if (!strcmp(name, "o2")) {
        sep_header->o2 = atoi(value);
      } else if (!strcmp(name, "o3")) {
        sep_header->o3 = atoi(value);
      } else if (!strcmp(name, "d1")) {
        sep_header->d1 = atof(value);
      } else if (!strcmp(name, "d2")) {
        sep_header->d2 = atof(value);
      } else if (!strcmp(name, "d3")) {
        sep_header->d3 = atof(value);
      } else if (!strcmp(name, "label1")) {
        strncpy(sep_header->label1, value, value_len+1);
      } else if (!strcmp(name, "label2")) {
        strncpy(sep_header->label2, value, value_len+1);
      } else if (!strcmp(name, "label3")) {
        strncpy(sep_header->label3, value, value_len+1);
      } else if (!strcmp(name, "unit1")) {
        strncpy(sep_header->unit1, value, value_len+1);
      } else if (!strcmp(name, "unit2")) {
        strncpy(sep_header->unit2, value, value_len+1);
      } else if (!strcmp(name, "unit3")) {
        strncpy(sep_header->unit3, value, value_len+1);
      } else if (!strcmp(name, "data_format")) {
        strncpy(sep_header->data_format, value, value_len+1);
      } else if (!strcmp(name, "in")) {
        strncpy(sep_header->in, value, value_len+1);
      }
    }

    len = 0;
  }
  if (line)
    free(line);
  fclose(fd);
  regfree(&re);
}

/*****************************************************************************/
void parse_sep_vs_binary(sep_header_t *vs_header, int *topo) {
  FILE *fd = fopen(vs_header->in, "r");
  if (fd == NULL)
    EXIT_ON_ERR("fopen");

  size_t buf_size = vs_header->n1 * vs_header->n2;
  float *buf = (float *) malloc(buf_size * sizeof(float));


  int count = 0;
  memset(topo, -1, vs_header->n1 * vs_header->n2 *sizeof(int));


  for (int k = 0; k < vs_header->n3; ++k) {
    size_t num_read = fread(buf, sizeof(float), buf_size, fd);
    if (num_read != buf_size)
      EXIT_ON_ERR("fread");
    for (int j = 0; j < vs_header->n2; ++j) {
      for (int i = 0; i < vs_header->n1; ++i) {
        if ((topo[i + j*vs_header->n1] == -1) &&
            (buf[i + j*vs_header->n1] != 0.0)) {
          topo[i + j*vs_header->n1] = k;
          count++;
        }
      }
    }
    if (count == vs_header->n1 * vs_header->n2)
      break;
  }
  free(buf);
  fclose(fd);
  fprintf(stderr, "Found : %d bathymetry points.\n", count);
}
