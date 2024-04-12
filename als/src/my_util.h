#pragma once


#include "header.h"


const double EPSILON = 1e-8;
const double DELAY_TOL = 1e-3;
const double AREA_TOL = 1e-3;


int ExecSystComm(const std::string & cmd);
bool IsPathExist(const std::string & _path);
void CreatePath(const std::string & _path);
void FixPath(std::string & _path);


static inline bool DoubleEqual     (const double a, const double b, double epsilon = EPSILON) {return fabs(a - b) < epsilon;}
static inline bool DoubleGreat     (const double a, const double b, double epsilon = EPSILON) {return a - b >= epsilon;}
static inline bool DoubleGreatEqual(const double a, const double b, double epsilon = EPSILON) {return a - b > -epsilon;}
static inline bool DoubleLess      (const double a, const double b, double epsilon = EPSILON) {return a - b <= -epsilon;}
static inline bool DoubleLessEqual (const double a, const double b, double epsilon = EPSILON) {return a - b < epsilon;}
static inline void FixAndCreatePath(std::string & _path) {FixPath(_path); CreatePath(_path);}


template <typename T>
void PrintVect(const std::vector <T> & vect, std::string && endWith = "", ll len = std::numeric_limits <ll>::max(), bool isRev = false) {
    std::cout << "(";
    ll size = std::min(len, static_cast <ll> (vect.size()));
    if (!isRev) {
        for (ll i = 0; i < size; ++i) {
            if (i == vect.size() - 1)
                std::cout << vect[i];
            else
                std::cout << vect[i] << ",";
        }
    }
    else {
        for (ll i = size - 1; i >= 0; --i) {
            if (i == 0)
                std::cout << vect[i];
            else
                std::cout << vect[i] << ",";
        }
    }
    std::cout << ")";
    std::cout << endWith;
}


template <typename T>
void PrintVect(const std::vector <T> && vect, std::string && endWith = "") {
    std::cout << "(";
    for (ll i = 0; i < vect.size(); ++i) {
        if (i == vect.size() - 1)
            std::cout << vect[i];
        else
            std::cout << vect[i] << ",";
    }
    std::cout << ")";
    std::cout << endWith;
}


template <typename T>
void PrintVect(const std::vector < std::vector <T> > & vect) {
    for (auto & line: vect)
        PrintVect <T> (line, "\n");
}


template <typename T>
concept INT_NUMB = (
    std::is_same_v<T, ull> ||
    std::is_same_v<T, ll>
);


template <INT_NUMB T>
inline void SetBit(T & x, const ll bit) {
    x |= (1llu << bit);
}


template <INT_NUMB T>
inline void ResetBit(T & x, const ll bit) {
    x &= ~(1llu << bit);
}


template <INT_NUMB T>
inline void ChangeBit(T & x, const ll bit, const bool val) {
    x |= (static_cast <ull> (val) << bit);
}


template <INT_NUMB T>
inline bool GetBit(const T x, const ll bit) {
    return static_cast <bool> ((x >> bit) & 1llu);
}


template <>
inline void
boost::to_block_range(const boost::dynamic_bitset<>& b, std::tuple<ll, ull&> param)
{
    ll iBlock = std::get<0>(param);
    std::get<1>(param) = b.m_bits[iBlock];
    return;
}


static inline ull GetBlockFromDynBitset(const boost::dynamic_bitset<>& b, ll iBlock) {
    ull res = 0;
    boost::to_block_range(b, std::make_tuple(iBlock, std::ref(res)));
    return res;
}