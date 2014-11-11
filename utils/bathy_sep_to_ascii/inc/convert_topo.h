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

#ifndef _CONVERT_TOPO_H_
#define _CONVERT_TOPO_H_

void convert_topo_sep_to_ascii(char *sep_name, char *ascii_name);

void read_topo_sep(sep_header_t *header, int *sep_topo);

void assign_ascii_topo_values(sep_header_t *header,
                              int *sep_topo, float *ascii_topo);

void write_ascii_topo(char *ascii_name, float *ascii_topo, int n);

#endif /* end of include guard: _CONVERT_TOPO_H_ */
