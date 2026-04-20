//
// Created by raj on 4/20/26.
//

#include "core/error.h"

#include <stdio.h>
#include <stdarg.h>
#include <string.h>

#define MAX_ERR_LEN 512 ///isse lamba kya hi hai andf in fact we woudl rather focus on the code
static _Thread_local char last_error_msg[MAX_ERR_LEN] = {0};

void clearerror() {
    last_error_msg[0] = '\0';
}

void setlasterror(const char *format, ...) {
    va_list args;
    va_start(args, format);

    vsnprintf(last_error_msg, MAX_ERR_LEN, format, args);

    va_end(args);

    fprintf(stderr, "[DEBUG ERROR]: %s\n", last_error_msg);  /// alkwayts dump the dfsmn error (will ifndef later)

}
const char *getlasterror() {
    return last_error_msg;
}