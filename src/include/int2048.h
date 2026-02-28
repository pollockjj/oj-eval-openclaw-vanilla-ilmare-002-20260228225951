#pragma once
#ifndef SJTU_BIGINTEGER
#define SJTU_BIGINTEGER

#include <complex>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

namespace sjtu {

class int2048 {
public:
    static const long long BASE = 100000000LL; // 10^8
    static const int BASEDIGITS = 8;
    std::vector<long long> d; // little-endian, d[0] is least significant
    bool neg;                  // true if negative

    void trim();
    bool isZero() const;
    int cmpAbs(const int2048 &o) const;
    void addAbs(const int2048 &o);
    void subAbs(const int2048 &o);

    static std::vector<long long> multiply_vec(const std::vector<long long> &a, const std::vector<long long> &b);
    static void divmod(const int2048 &a, const int2048 &b, int2048 &q, int2048 &r);

public:
    int2048();
    int2048(long long);
    int2048(const std::string &);
    int2048(const int2048 &);

    void read(const std::string &);
    void print();

    int2048 &add(const int2048 &);
    friend int2048 add(int2048, const int2048 &);

    int2048 &minus(const int2048 &);
    friend int2048 minus(int2048, const int2048 &);

    int2048 operator+() const;
    int2048 operator-() const;

    int2048 &operator=(const int2048 &);

    int2048 &operator+=(const int2048 &);
    friend int2048 operator+(int2048, const int2048 &);

    int2048 &operator-=(const int2048 &);
    friend int2048 operator-(int2048, const int2048 &);

    int2048 &operator*=(const int2048 &);
    friend int2048 operator*(int2048, const int2048 &);

    int2048 &operator/=(const int2048 &);
    friend int2048 operator/(int2048, const int2048 &);

    int2048 &operator%=(const int2048 &);
    friend int2048 operator%(int2048, const int2048 &);

    friend std::istream &operator>>(std::istream &, int2048 &);
    friend std::ostream &operator<<(std::ostream &, const int2048 &);

    friend bool operator==(const int2048 &, const int2048 &);
    friend bool operator!=(const int2048 &, const int2048 &);
    friend bool operator<(const int2048 &, const int2048 &);
    friend bool operator>(const int2048 &, const int2048 &);
    friend bool operator<=(const int2048 &, const int2048 &);
    friend bool operator>=(const int2048 &, const int2048 &);
};

} // namespace sjtu

#endif
