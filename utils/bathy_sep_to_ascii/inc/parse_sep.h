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
#include "sep_header.h"

#ifndef _SEP_PARSER_H_
#define _SEP_PARSER_H_

/*****************************************************************************/

/**
  * \brief Fill a sep_header structure from a sep_header
  *
  * \param vs_name Name of the sep header
  * \param sep_header Structure to fill
  */
void parse_sep_vs_header(char *vs_name, sep_header_t *sep_header);

/**
  * \brief Get bathimetry location, where vs == 0.
  *
  * \param vs_header SEP header indicating dimensions and file location.
  * \param topo Array to store index where liquid / solid interface is found.
  */
void parse_sep_vs_binary(sep_header_t *vs_header, int *topo);

/*****************************************************************************/

#endif /* end of include guard: _SEP_PARSER_H_ */
