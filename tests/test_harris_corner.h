#if !defined(_TEST_HARRIS_CORNER_H_)
#define _TEST_HARRIS_CORNER_H_

#include "test_config.h"

#if TEST_AVAILABLE_ALL

#include "mxm/cv_corner.h"
#include "mxm/random.h"
#endif

namespace mxm{

inline void testHarrisCorner()
{
#if TEST_AVAILABLE_ALL
    Mat img = random::uniform<float>({10,10});
    auto ret = harrisCornernessMap(img, 5);
#endif
}

} //scope of namespace mxm
#endif // _TEST_HARRIS_CORNER_H_
