#include <iostream>
#include "test_aabb.h"
// #include "test_rigid_body.h"
#include "test_linalg.h"
#include "test_ray.h"
#include "test_random.h"
// #include "test_material.h"
#include "test_rotation.h"
#include "test_bvh.h"
#include "test_interp.h"
#include "test_rigid_transform.h"
#include "test_bazier.h"
#include "test_camera.h"
// #include "test_pixel.h"
// #include "test_image_processing.h"
#include "test_bsp_tree.h"


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
  testBspTree();
  std::cout << "done" << std::endl;
  return 0;
}
