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

#ifndef _SEP_HEADER_H_
#define _SEP_HEADER_H_

/*****************************************************************************/

/**
  * brief SEP header informations
  *
  * 3D files at most.
  */
typedef struct {
  int n1, n2, n3;
  int o1, o2, o3;
  float d1, d2, d3;
  char label1[512], label2[512], label3[512];
  char unit1[512], unit2[512], unit3[512];
  char data_format[512];
  char in[512];
} sep_header_t;

/**
  * \brief Display the values in a sep header structure.
  *
  * \param h Header information to display
  */
void show_header(sep_header_t h);

/*****************************************************************************/

#endif /* end of include guard: _SEP_HEADER_H_ */
