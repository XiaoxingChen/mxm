#if !defined(_HAMMING_SPACE_H_)
#define _HAMMING_SPACE_H_

#include <bitset>
#include <string>
#include <sstream>
#include <vector>
#include <map>

namespace mxm
{

template<size_t N>
std::string to_string(const std::bitset<N>& val)
{
    std::stringstream stream;
    stream << std::setfill ('0') << std::uppercase;
    stream << std::hex;
    stream << val.to_ullong();
    // for(size_t i = 0; i < (N / 8); i++)
    // {
    //     size_t pos = K - 1 - i;
    //     stream << std::setw(2) << ((val.to_ullong() >> (pos * 8)) & 0xff);
    // }
    return stream.str();
}

template<size_t N>
std::string to_string(const std::bitset<N>& val, size_t prec)
{
    return to_string<N>(val);
}

#if 1
template<size_t N>
std::multimap<size_t, size_t> nearestNeighborSearch(const std::vector<std::bitset<N>>& pool, const std::bitset<N>& val, size_t k)
{
    std::multimap<size_t, size_t> ret;
    for(size_t i = 0; i < pool.size(); i++)
    {
        size_t dist = (pool.at(i) ^ val).count();
        if(dist < std::prev(ret.end())->first)
        {
            ret.insert({dist, i});
        }

        if(ret.size() > k)
        {
            ret.erase(std::prev(ret.end()));
        }
    }
    return ret;
}
#endif

template<size_t N>
std::multimap<size_t, size_t>
radiusSearch(const std::vector<std::bitset<N>>& pool, const std::bitset<N>& val, size_t r)
{
    std::multimap<size_t, size_t> ret;
    for(size_t i = 0; i < pool.size(); i++)
    {
        size_t dist = (pool.at(i) ^ val).count();
        if(dist <= r)
        {
            ret.insert({dist, i});
        }
    }
    return ret;
}


} // namespace mxm

#endif // _HAMMING_SPACE_H_
