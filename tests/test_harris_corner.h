#if !defined(_TEST_HARRIS_CORNER_H_)
#define _TEST_HARRIS_CORNER_H_

#include "mxm/cv_corner.h"
#include "mxm/random.h"

using namespace mxm;

inline void testHarrisCorner()
{
    Mat img = random::uniform<float>({10,10});
    auto ret = harrisCornernessMap(img, 5);
}


#endif // _TEST_HARRIS_CORNER_H_
