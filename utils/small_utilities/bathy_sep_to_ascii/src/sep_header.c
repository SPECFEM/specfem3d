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
#include "sep_header.h"

/*****************************************************************************/
void show_header(sep_header_t h) {
  fprintf(stderr, "n1: %d / n2: %d / n3: %d\n", h.n1, h.n2, h.n3);
  fprintf(stderr, "o1: %d / o2: %d / o3: %d\n", h.o1, h.o2, h.o3);
  fprintf(stderr, "d1: %f / d2: %f / d3: %f\n", h.d1, h.d2, h.d3);
  fprintf(stderr, "label1: %s / label2: %s / label3: %s\n", h.label1, h.label2, h.label3);
  fprintf(stderr, "unit1: %s / unit2: %s / unit3: %s\n", h.unit1, h.unit2, h.unit3);
  fprintf(stderr, "data_format: %s\n", h.data_format);
  fprintf(stderr, "in: %s\n", h.in);
}
