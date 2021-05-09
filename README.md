# mxm

A header-only/compiled C++ numerical compute library.

<img src="https://render.githubusercontent.com/render/math?math=A_{m\times m}">

# Features

## 1. Matrix
```cpp
#include "mxm/linalg.h"
using namespace mxm;
```
- Initialize a matrix

```cpp
Mat mat_a = Mat({3,3}, {
    1,1,1,
    2,2,2,
    3,3,3})
```

- Matrix block assignment

```cpp
Mat transform = Mat({4,4});
transform(3,3) = 1;
transform(Block({0,3},{0,3})) = Mat::Identity(3);
```

- Random matrix

```cpp
Mat mat_a = random::uniform({5,5});
```

- Common operations: matrix multiplication, transpose, inversion

```cpp
Mat mat_a = random::uniform({3,3});
Mat mat_b = Mat::Identity({3,3});
Mat mat_c = mat_a.matmul(mat_b);
mat_a = mat_a.T() // O(1) time complexity
```

```cpp
Mat mat_a = Mat::ones({3,3});
mat_a(Col(2)) += Vec({2,3,4});
Mat inv_b = mat_a.inv();

```

- QR decomposition and linear equation

```cpp
Mat mat_a = random::uniform({5,5});
Vec vec_b = random::uniform({5,1});
auto qr = qr::decomposeByRotation(mat_a);
Vec vec_x = qr::solve(mat_a, vec_b);
```

## 2. N-Dimensional Rigid Body Transform

```cpp
#include "mxm/rigid_transform.h"
using namespace mxm;
```

- 2D Rotation
```cpp
Rotation r1 = Rotation::fromAngle(0.1);
Rotation r2 = Rotation::fromAngle(0.5);
Vec v1 = (r1 * r2).apply(Vec({1,2}));
```

- 3D Rotation

```cpp
Rotation r1 = Rotation::fromAxisAngle(Vec({0,0,1}), 0.5);
Rotation r2 = Rotation::fromAxisAngle(Vec({0,1,0}), 0.3);
Vec v1 = (r1 * r2).apply(Vec({1,2,3}));
```

- 4D Rotation

```cpp
// Rotation::fromPlaneAngle() is available for any dimensional rotation.
Rotation r1 = Rotation::fromPlaneAngle(Vec({1,0,0,0}), Vec({0,1,0,0}), 0.5);
Vec v1 = r1.apply(Vec({1,2,3,4}));
```

- 3D Rigid Body Transform

```cpp
RigidTrans tf = RigidTrans::Identity(3);
Vec v1 = tf.apply(Vec({1,0,0}));
```

## 3. Space Indexing

```cpp
#include "mxm/spatial_bvh.h"
using namespace mxm;
```

- BVH Tree and Ray-primitive Intersection

```cpp
Mat vertex_buffer;
Matrix<size_t> vertex_buffer;
createTriangleBand(&vertex_buffer, &index_buffer);
bvh::PrimitiveMeshTree tree(vertex_buffer, index_buffer);
tree.build(4, false);

Ray ray({0,0,0},{1,0,0});
auto records = tree.hit(ray, bvh::eMultiHit);
```

- Radius Search and K Nearest Neighbor Search

```cpp
size_t dim = 2;
std::shared_ptr<Mat> pts(new Mat(random::uniform<FloatType>({dim, 100})));

bvh::PointCloudTree tree(pts);
tree.build(4, false);

Vec target_pt({.5, .5});
FloatType radius = 0.5;

auto radius_search_result = tree.radiusSearch(target_pt, radius);

size_t k = 5;
std::multimap<FloatType, size_t> knn_result = tree.nearestNeighborSearch(target_pt, k);
```



# Projects using mxm
1. [ray_tracing_4d](https://github.com/XiaoxingChen/ray_tracing_4d)