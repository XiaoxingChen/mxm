#include <iostream>
#include "test_linalg.h"
#include "test_aabb.h"
// #include "test_rigid_body.h"
#include "test_ray.h"
#include "test_random.h"
// #include "test_material.h"
#include "test_rotation.h"
#include "test_bvh.h"
#include "test_interp.h"
#include "test_rigid_transform.h"
#include "test_bazier.h"
#include "test_camera.h"
#include "test_cv_basic.h"
// #include "test_image_processing.h"
#include "test_bsp_tree.h"
#include "test_graph.h"
#include "test_grid_map.h"
#include "test_harris_corner.h"
#include "test_optical_flow.h"
#include "test_hamming_space.h"
#include "test_coordinate_system.h"
#include "test_lie_special_orthogonal.h"
#include "test_lie_special_unitary.h"

int main(int argc, char const *argv[])
{
  testAABB();
  testRay();
  testLinearAlgebra();
  // testPixel();
  testRotation();
  testInterpolation();
  // testRigidBody();
  // testPrimitiveGeometry();
  testRandom();
  // testMaterial();
  testBvh();
  testRigidTransform();
  testBazier();
  testCamera();
  // testImageProcessing();
  testCvBasic();
  testBspTree();
  testGraph();
  testGridMap();
  testHarrisCorner();
  testOpticalFlow01();
  testHammingSpace();
  testCoordinateSystem();
  testLieSpecialOrthogonal();
  testLieSpecialUnitary();

  std::cout << "done" << std::endl;
  return 0;
}
