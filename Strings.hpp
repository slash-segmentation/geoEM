#pragma once

#include <string>
#include <algorithm>

#include <cctype>
#include <cstring>

// C++ std::string functions
inline std::string tolower(std::string s) { std::transform(s.begin(), s.end(), s.begin(), (int(*)(int))std::tolower); return s; }
inline void tolower_inplace(std::string& s) { std::transform(s.begin(), s.end(), s.begin(), (int(*)(int))std::tolower); }
inline std::string trim(const std::string& s)
{
    size_t i = 0, j = s.size() - 1;
    while (i <= j && std::isspace(s[i])) { ++i; }
    while (i <= j && std::isspace(s[j])) { --j; }
    return s.substr(i, j+1-i);
}
inline void trim_inplace(std::string& s)
{
    size_t i = 0, j = s.size() - 1;
    while (i <= j && std::isspace(s[i])) { ++i; }
    while (i <= j && std::isspace(s[j])) { --j; }
    s.erase(j+1).erase(0, i);
}
inline ssize_t read_int(const std::string& s)
{
    if (s.size() == 0) { throw std::invalid_argument("Error: invalid integer"); }
    size_t i = 0;
    bool negative = s[0] == '-';
    if (negative || s[0] == '+') { if (s.size() == 1) { throw std::invalid_argument("Error: invalid integer"); } i = 1; }
    ssize_t _x = 0;
    do
    {
        if (!std::isdigit(s[i])) { throw std::invalid_argument("Error: invalid integer"); }
        int64_t temp = ((int64_t)_x) * 10 + (s[i++] - '0');
        if (temp > INT_MAX) { throw std::invalid_argument("Error: invalid integer"); } // overflowed
        _x = (int)temp;
    }
    while (i < s.size());
    return negative ? -_x : _x;
}

// C-style string functions
inline char *trim(char *s) { while (std::isspace(*s)) { ++s; } size_t l = strlen(s); while (l && std::isspace(s[l - 1])) { s[--l] = 0; } return s; }
inline char *ltrim(char *s) { while (std::isspace(*s)) { ++s; } return s; }
inline const char *ltrim(const char *s) { while (std::isspace(*s)) { ++s; } return s; }
inline char *strchr(char *s, const char *cs) { while (*s) { const char *c = cs; while (*c) { if (*s == *c) { return s; } ++c; } ++s; } return nullptr; }
inline const char *strchr(const char *s, const char *cs) { while (*s) { const char *c = cs; while (*c) { if (*s == *c) { return s; } ++c; } ++s; } return nullptr; }
inline bool streq(const char *a, const char *b) { return std::strcmp(a, b) == 0; }
inline bool strieq(const char *a, const char *b)
{
    while (*a && *b && std::tolower(*a) == std::tolower(*b)) { ++a; ++b; }
    return *a == 0 && *b == 0;
}
inline bool str_read_int(const char *s, char **endptr, ssize_t *x)
{
    bool negative = *s == '-';
    if (negative || *s == '+') { ++s; }
    if (!std::isdigit(*s)) { return false; }
    ssize_t _x = 0;
    do
    {
        int64_t temp = ((int64_t)_x) * 10 + (*s++ - '0');
        if (temp > INT_MAX) { return false; } // overflowed
        _x = (int)temp;
    }
    while (std::isdigit(*s));
    *x = (negative ? -_x : _x);
    *endptr = (char*)s;
    return true;
}
inline bool str_skip_int(const char *s, char **endptr)
{
    if (*s == '-' || *s == '+') { ++s; }
    if (!std::isdigit(*s)) { return false; }
    while (std::isdigit(*++s));
    *endptr = (char*)s;
    return true;
}
