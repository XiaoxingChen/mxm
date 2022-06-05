# mxm

A header-only/compiled C++ numerical compute library.

<!-- <img src="https://render.githubusercontent.com/render/math?math=A_{m\times m}"> -->

$$\Huge A_{m \times m}$$

# Installation

## Linux/MacOS

```sh
git clone https://github.com/XiaoxingChen/mxm
cd mxm
./build.py --test
```
## Windows
```sh
git clone https://github.com/XiaoxingChen/mxm
cd mxm
python build.py --test
```
# Features

### 1. Dense Matrix Operation

- N-Dimensional Matrix/Vector Operation [<a href="./inc/mxm/linalg_mat.h">mxm/linalg_mat.h</a>]
    - Real
    - Complex [<a href="./inc/mxm/linalg_complex.h">mxm/linalg_complex.h</a>]
    - Quaternion [<a href="./inc/mxm/linalg_complex.h">mxm/linalg_complex.h</a>]
    - Dual Number [<a href="./inc/mxm/linalg_dual_number.h">mxm/linalg_dual_number.h</a>]
- QR Decomposition/solver (Real, Complex) [<a href="./inc/mxm/linalg_solve.h">mxm/linalg_solve.h</a>]
- Eigenvalue Decomposition (Real, by shifted QR Iteration).
- Singular Value Decomposition (Real).

### 2. Geometry

- 2-D, 3-D, N-D Rotation
    - Logarithmic/Exponential Map for $SO(n)$ and $\mathfrak{so}(n)$ [<a href="./inc/mxm/lie_special_orthogonal.h">mxm/lie_special_orthogonal.h</a>]
    - Conversion between Axis/Plane-Angle, Rotation Matrix, Quaternion [<a href="./inc/mxm/full_dimensional_rotation.h">mxm/full_dimensional_rotation.h</a>]
    - Spherical Linear Interpolation
- N-D Ray [<a href="./inc/mxm/geometry_ray.h">mxm/geometry_ray.h</a>]
- N-D Pinhole Camera (3D Radial Tangential Distortion) [<a href="./inc/mxm/model_camera.h">mxm/model_camera.h</a>]
- N-D Simplex(line segment, triangle, tetrahedron ...) [<a href="./inc/mxm/geometry_primitive.h">mxm/geometry_primitive.h</a>]
- N-D Rigidbody Transform [<a href="./inc/mxm/rigid_transform.h">mxm/rigid_transform.h</a>]
- N-D Affine Transform [<a href="./inc/mxm/transform_affine.h">mxm/transform_affine.h</a>]
- N-D BVH Tree [<a href="./inc/mxm/spatial_bvh.h">mxm/spatial_bvh.h</a>]
    - Ray Closest-Hit/Any-Hit/Multi-Hit
    - Radius Search
    - Nearest Neighbor Search
- N-D Metric Tree (for any metric space) [<a href="./inc/mxm/spatial_metric_tree.h">mxm/spatial_metric_tree.h</a>]

### 3. None-linear Optimization

- Gauss-Newton Method [<a href="./inc/mxm/optimize.h">mxm/optimize.h</a>]
- Auto Derivative(based on Dual Number) :star2: :star2:
- Sparse Jacobian Acceleration for SfM [<a href="./inc/mxm/optimize_sfm_gn.h">mxm/optimize_sfm_gn.h</a>]

### 4. Graph

- Directed/Undirected Weighted/Unweighted Graph [<a href="./inc/mxm/graph_base.h">mxm/graph_base.h</a>]
- Shortest Path [<a href="./inc/mxm/graph_dijkstra.h">mxm/graph_dijkstra.h</a>]
    - Dijkstra
    - Bellman Ford
- Flow Network [<a href="./inc/mxm/graph_flow.h">mxm/graph_flow.h</a>]
    - Fulkerson Ford Max Flow

### 5. Toy Demos
- Optical Flow [<a href="./inc/mxm/cv_optical_flow.h">mxm/cv_optical_flow.h</a>]
- Pinhole Camera Calibration [<a href="./inc/mxm/cv_calibration.h">mxm/cv_calibration.h</a>]


# Code Sample
- QR decomposition and linear equation

```cpp
#include "mxm/linalg.h"
using namespace mxm;

Matrix<float> mat_a = random::uniform({5,5});
Vector<float> vec_b = random::uniform({5,1});
auto qr = qr::decomposeByRotation(mat_a);
auto vec_x = qr::solve(mat_a, vec_b);
```

- 4D Rotation

```cpp
// Rotation::fromPlaneAngle() is available for any dimensional rotation.
Rotation<double> r1 = Rotation::fromPlaneAngle({1,0,0,0}, {0,1,0,0}, 0.5);
auto v1 = r1.apply({1,2,3,4});
```

- BVH Tree and Ray-primitive Intersection

```cpp
#include "mxm/spatial_bvh.h"
using namespace mxm;

Mat vertex_buffer;
Matrix<size_t> vertex_buffer;
createTriangleBand(&vertex_buffer, &index_buffer);
bvh::PrimitiveMeshTree tree(vertex_buffer, index_buffer);
tree.build(4, false);

Ray ray({0,0,0},{1,0,0});
auto records = tree.hit(ray, bvh::eMultiHit);
```

For more code samples, see [tests/test_main.cpp](./tests/test_main.cpp).

# Related Blogs and Links
1. A 4D CPU Ray Tracing Renderer based on `mxm`: [ray_tracing_4d](https://github.com/XiaoxingChen/ray_tracing_4d)
2. Blog: [Dual Number and Auto Derivative of Special Orthogonal Group](https://xiaoxingchen.github.io/2022/03/20/dual_number_auto_derivative_on_so3/)
3. Blog: [Geometries for N-Dimensional Ray Tracing](https://xiaoxingchen.github.io/2021/01/07/ray_tracing_4d_01/)
4. Zhihu: [How to implement Auto Differentiation in C++?](https://www.zhihu.com/question/48356514/answer/2446699680)(Chinese)
