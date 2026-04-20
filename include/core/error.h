//
// Created by raj on 4/20/26.
//

#ifndef QISC_ERROR_H
#define QISC_ERROR_H

void clearerror(void);
void setlasterror(const char *format, ...);
const char *getlasterror(void);
#endif //QISC_ERROR_H
