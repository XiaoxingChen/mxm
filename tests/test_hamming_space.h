#if !defined(_TEST_HAMMING_SPACE_H_)
#define _TEST_HAMMING_SPACE_H_

#include <iostream>

#include "mxm/hamming_space.h"
#include "test_config.h"

namespace mxm
{
#if TEST_HAMMING_SPACE
void testHammingSpace()
{

{
    std::vector<std::bitset<32>> pool{0xaabbccdd, 0x01010202, 0xbc030405};
    auto result = radiusSearch(pool, std::bitset<32>(0xabbbcddd), 5);
    if(result.begin()->second != 0)
    {
        std::cout << mxm::to_string(pool[0]) << std::endl;
        throw std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__));
    }


}

}
#else
void testHammingSpace() {}
#endif


} // namespace mxm

#endif // _TEST_HAMMING_SPACE_H_
