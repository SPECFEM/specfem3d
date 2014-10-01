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
#include "err.h"

#define EXIT_ON_ERR(msg) \
    err(1, "%s(%d) -- %s: %s", __FILE__, __LINE__, __func__, msg);


/*****************************************************************************/
void parse_sep_header(char *header_name, 
                      int *n1, int *n2, int *n3, 
                      float *o1, float *o2, float *o3, 
                      float *d1, float *d2, float *d3, 
                      char *in) {
  char *line = NULL;
  size_t len = 0;

  FILE *fd = fopen(header_name, "r");
  if (fd == NULL) 
    EXIT_ON_ERR(header_name);

  /* Build regex mathing "name=value" with optional spaces*/
  regex_t re;
  if(regcomp(&re, "^[ \t]*\\(.*\\)[ \t]*=[ \t]*\\(.*\\)[ \t]*", 
             REG_ICASE | REG_NEWLINE))
    EXIT_ON_ERR("regcomp");

  /* Apply the regex to every line of the sep header */
  while(getline(&line, &len, fd) != -1) {

    char name[512], value[512];
    size_t nmatch = 3;
    regmatch_t pmatch[3];
    if(REG_NOMATCH == regexec(&re, line, nmatch, pmatch, 0)) {
      continue; /* These files are usually quite messy... */
    } else {
      size_t name_len = pmatch[1].rm_eo - pmatch[1].rm_so;
      strncpy(name, &line[pmatch[1].rm_so], name_len);
      name[name_len] = '\0';

      size_t value_len = pmatch[2].rm_eo - pmatch[2].rm_so;
      strncpy(value, &line[pmatch[2].rm_so], pmatch[2].rm_eo - pmatch[2].rm_so);
      value[value_len] = '\0';

      if (!strcmp(name, "n1")) {
        *n1 = atoi(value);
      } else if (!strcmp(name, "n2")) {
        *n2 = atoi(value);
      } else if (!strcmp(name, "n3")) {
        *n3 = atoi(value);
      } else if (!strcmp(name, "o1")) {
        *o1 = atof(value);
      } else if (!strcmp(name, "o2")) {
        *o2 = atof(value);
      } else if (!strcmp(name, "o3")) {
        *o3 = atof(value);
      } else if (!strcmp(name, "d1")) {
        *d1 = atof(value);
      } else if (!strcmp(name, "d2")) {
        *d2 = atof(value);
      } else if (!strcmp(name, "d3")) {
        *d3 = atof(value);
      /*} else if (!strcmp(name, "label1")) {
        strncpy(*label1, value, value_len+1);
      } else if (!strcmp(name, "label2")) {
        strncpy(*label2, value, value_len+1);
      } else if (!strcmp(name, "label3")) {
        strncpy(*label3, value, value_len+1);
      } else if (!strcmp(name, "unit1")) {
        strncpy(*unit1, value, value_len+1);
      } else if (!strcmp(name, "unit2")) {
        strncpy(*unit2, value, value_len+1);
      } else if (!strcmp(name, "unit3")) {
        strncpy(*unit3, value, value_len+1);
      } else if (!strcmp(name, "data_format")) {
        strncpy(*data_format, value, value_len+1);*/
      } else if (!strcmp(name, "in")) {
        strncpy(in, value, value_len+1);
      }
    }

    if (line)
      free(line);
    len = 0;
  }
  fclose(fd);
  regfree(&re);
}


/*****************************************************************************/
/*** Aliases for Fortran bindings ***/

void parse_sep_header_(char *header_name, 
                       int *n1, int *n2, int *n3, 
                       int *o1, int *o2, int *o3, 
                       float *f1, float *f2, float *f3, 
                       char *in)
    __attribute__((alias("parse_sep_header")));
/*****************************************************************************/
