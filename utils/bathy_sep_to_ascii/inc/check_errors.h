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

#ifndef _CHECK_ERRORS_H_
#define _CHECK_ERRORS_H_

/**
 * \brief Macro to be called in case a function return an error code,
 *        print a message on stderr and exit.
 */
#ifdef HAVE_ERR
#include "err.h"
#define EXIT_ON_ERR(msg) \
    err(1, "%s(%d) -- %s: %s", __FILE__, __LINE__, __func__, msg);
#else
#define EXIT_ON_ERR(msg) \
    printf("%s(%d) -- %s: %s", __FILE__, __LINE__, __func__, msg); exit(1);
#endif

#endif /* end of include guard: _CHECK_ERRORS_H_ */
