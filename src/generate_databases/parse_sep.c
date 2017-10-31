/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  3 . 0
 !               ---------------------------------------
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

/*****************************************************************************/
#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>
#ifdef HAVE_ERR
#include "err.h"
#define EXIT_ON_ERR(msg) \
    err(1, "%s(%d) -- %s: %s", __FILE__, __LINE__, __func__, msg);
#else
#define EXIT_ON_ERR(msg) \
    printf("%s(%d) -- %s: %s", __FILE__, __LINE__, __func__, msg); exit(1);
#endif

/*****************************************************************************/
/**
 * Parse a sep header file and get information of interest
 *
 * \param header_name Name of the sep header to parse
 * \param n1 Number of sample in X direction
 * \param n2 Number of sample in Y direction
 * \param n3 Number of sample in Z direction
 * \param o1 Offset (m) in X direction
 * \param o2 Offset (m) in Y direction
 * \param o3 Offset (m) in Z direction
 * \param d1 Spatial increment (m) in X direction
 * \param d2 Spatial increment (m) in Y direction
 * \param d3 Spatial increment (m) in Z direction
 * \param in Name and location of the SEP binary file
 */
void
FC_FUNC_(parse_sep_header, PARSE_SEP_HEADER)
    (char *header_name,
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
  while (getline(&line, &len, fd) != -1) {

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
        /* We do not have use for other fields. */
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
  }

  // buffer will be reallocated by getline
  // should be freed even if getline is not successful
  free(line);

  fclose(fd);
  regfree(&re);
}
