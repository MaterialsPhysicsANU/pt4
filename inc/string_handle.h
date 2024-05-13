#ifndef STRING_HANDLE_H
#define STRING_HANDLE_H

#ifdef __cplusplus
#include<string>

typedef std::string std__string;

extern "C" {
#else
typedef struct std__string std__string;
#endif

typedef std__string * string_handle;
string_handle create_string(char*);
string_handle create_quoted_string(char*);
void push_back_quoted_string(string_handle p, char* s);
void    free_string(string_handle);
const char*    c_str(string_handle);

#ifdef __cplusplus
}
#endif

#endif