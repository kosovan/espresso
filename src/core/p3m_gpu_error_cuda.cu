#include <cuda.h>
#include <stdio.h>
#include <thrust/device_vector.h>
#include <thrust/reduce.h>

#include "cuda_utils.hpp"
#include "p3m_gpu_error.hpp"

#include "p3m_gpu_common.hpp"

#if defined(OMPI_MPI_H) || defined(_MPI_H)
#error CU-file includes mpi.h! This should not happen!
#endif

/** @TODO: Extend to hight order. This comes from some 1/sin expansion in
 * Hockney/Eastwood */

template <int cao>
__device__ static double p3m_analytic_cotangent_sum(int n, double mesh_i) {
  const double c = sqr(cos(PI * mesh_i * n));

  switch (cao) {
  case 1:
    return 1;
  case 2:
    return (1.0 + c * 2.0) / 3.0;
  case 3:
    return (2.0 + c * (11.0 + c * 2.0)) / 15.0;
  case 4:
    return (17.0 + c * (180.0 + c * (114.0 + c * 4.0))) / 315.0;
  case 5:
    return (62.0 + c * (1072.0 + c * (1452.0 + c * (247.0 + c * 2.0)))) /
           2835.0;
  case 6:
    return (1382.0 +
            c * (35396.0 +
                 c * (83021.0 + c * (34096.0 + c * (2026.0 + c * 4.0))))) /
           155925.0;
  case 7:
    return (21844.0 +
            c * (776661.0 +
                 c * (2801040.0 +
                      c * (2123860.0 +
                           c * (349500.0 + c * (8166.0 + c * 4.0)))))) /
           6081075.0;
  }

  return 0.0;
}

template <int cao>
__global__ void p3m_k_space_error_gpu_kernel_ik(int3 mesh, double3 meshi,
                                                double alpha_L, double *he_q) {
  const int nx = -mesh.x / 2 + blockDim.x * blockIdx.x + threadIdx.x;
  const int ny = -mesh.y / 2 + blockDim.y * blockIdx.y + threadIdx.y;
  const int nz = -mesh.z / 2 + blockDim.z * blockIdx.z + threadIdx.z;

  if ((nx >= mesh.x / 2) || (ny >= mesh.y / 2) || (nz >= mesh.z / 2))
    return;

  const int lind = (nx + mesh.x / 2) * mesh.y * mesh.z +
                   (ny + mesh.y / 2) * mesh.z + (nz + mesh.z / 2);

  if ((nx != 0) || (ny != 0) || (nz != 0)) {
    const double alpha_L_i = 1. / alpha_L;
    const double n2 = sqr(nx) + sqr(ny) + sqr(nz);
    const double cs = p3m_analytic_cotangent_sum<cao>(nz, meshi.z) *
                      p3m_analytic_cotangent_sum<cao>(nx, meshi.x) *
                      p3m_analytic_cotangent_sum<cao>(ny, meshi.y);
    const double ex = exp(-(PI * alpha_L_i) * (PI * alpha_L_i) * n2);
    const double ex2 = sqr(ex);
    const double U2 = int_pow<2 * cao>(
        csinc(meshi.x * nx) * csinc(meshi.y * ny) * csinc(meshi.z * nz));
    auto const alias1 = ex2 / n2;
    auto const d = alias1 - sqr(U2 * ex / cs) / n2;

    if (d > 0 && (d / alias1 > ROUND_ERROR_PREC))
      he_q[lind] = d;
  } else {
    he_q[lind] = 0;
  }
}

__global__ void p3m_k_space_error_gpu_kernel_ad(int3 mesh, double3 meshi,
                                                int cao, double alpha_L,
                                                double *he_q) {
  int nx = -mesh.x / 2 + blockDim.x * blockIdx.x + threadIdx.x;
  int ny = -mesh.y / 2 + blockDim.y * blockIdx.y + threadIdx.y;
  int nz = -mesh.z / 2 + blockDim.z * blockIdx.z + threadIdx.z;

  if ((nx >= mesh.x / 2) || (ny >= mesh.y / 2) || (nz >= mesh.z / 2))
    return;

  int lind = ((nx + mesh.x / 2) * mesh.y * mesh.z + (ny + mesh.y / 2) * mesh.z +
              (nz + mesh.z / 2));

  double alpha_L_i = 1. / alpha_L;
  double n2;
  double U2, ex, ex2;
  int nmx, nmy, nmz;
  double alias1, alias2, alias3, alias4;

  alias1 = alias2 = alias3 = alias4 = 0;

  if ((nx != 0) || (ny != 0) || (nz != 0)) {
    for (int mx = -1; mx <= 1; mx++) {
      nmx = nx + mx * mesh.x;
      for (int my = -1; my <= 1; my++) {
        nmy = ny + my * mesh.y;
        for (int mz = -1; mz <= 1; mz++) {
          nmz = nz + mz * mesh.z;

          n2 = sqr(nmx) + sqr(nmy) + sqr(nmz);

          ex = exp(-(PI * alpha_L_i) * (PI * alpha_L_i) * n2);

          ex2 = sqr(ex);

          U2 = pow((double)csinc(meshi.x * nmx) * csinc(meshi.y * nmy) *
                       csinc(meshi.z * nmz),
                   2.0 * cao);

          alias1 += ex2 / n2;
          alias2 += U2 * ex;
          alias3 += U2 * n2;
          alias4 += U2;
        }
      }
    }

    if ((alias3 == 0.0) || (alias4 == 0.0))
      he_q[lind] = 0;
    else
      he_q[lind] = alias1 - (alias2 * alias2) / (alias3 * alias4);

  } else {
    he_q[lind] = 0;
  }
}

__global__ void p3m_k_space_error_gpu_kernel_ik_i(int3 mesh, double3 meshi,
                                                  int cao, double alpha_L,
                                                  double *he_q) {

  int nx = -mesh.x / 2 + blockDim.x * blockIdx.x + threadIdx.x;
  int ny = -mesh.y / 2 + blockDim.y * blockIdx.y + threadIdx.y;
  int nz = -mesh.z / 2 + blockDim.z * blockIdx.z + threadIdx.z;

  if ((nx >= mesh.x / 2) || (ny >= mesh.y / 2) || (nz >= mesh.z / 2))
    return;

  int lind = ((nx + mesh.x / 2) * mesh.y * mesh.z + (ny + mesh.y / 2) * mesh.z +
              (nz + mesh.z / 2));

  double alpha_L_i = 1. / alpha_L;
  double n2;
  double U2, ex, ex2;
  int nmx, nmy, nmz;
  double alias1, alias2, alias3, alias4;

  alias1 = alias2 = alias3 = alias4 = 0;

  if ((nx != 0) || (ny != 0) || (nz != 0)) {
    for (int mx = -1; mx <= 1; mx++) {
      nmx = nx + mx * mesh.x;
      for (int my = -1; my <= 1; my++) {
        nmy = ny + my * mesh.y;
        for (int mz = -1; mz <= 1; mz++) {
          nmz = nz + mz * mesh.z;

          n2 = sqr(nmx) + sqr(nmy) + sqr(nmz);

          ex = exp(-(PI * alpha_L_i) * (PI * alpha_L_i) * n2);

          ex2 = sqr(ex);

          U2 = pow((double)csinc(meshi.x * nmx) * csinc(meshi.y * nmy) *
                       csinc(meshi.z * nmz),
                   2.0 * cao);

          alias1 += ex2 / n2;
          alias2 += U2 * ex * (nx * nmx + ny * nmy + nz * nmz) / n2;
          alias3 += U2;

          if (((mx + my + mz) % 2) == 0) { // consider only even terms!
            alias4 += U2;
          } else {
            alias4 -= U2;
          }
        }
      }
    }

    he_q[lind] = alias1 -
                 (alias2 * alias2) / (0.5 * (nx * nx + ny * ny + nz * nz) *
                                      (alias3 * alias3 + alias4 * alias4));

  } else {
    he_q[lind] = 0;
  }
}

__global__ void p3m_k_space_error_gpu_kernel_ad_i(int3 mesh, double3 meshi,
                                                  int cao, double alpha_L,
                                                  double *he_q) {

  int nx = -mesh.x / 2 + blockDim.x * blockIdx.x + threadIdx.x;
  int ny = -mesh.y / 2 + blockDim.y * blockIdx.y + threadIdx.y;
  int nz = -mesh.z / 2 + blockDim.z * blockIdx.z + threadIdx.z;

  if ((nx >= mesh.x / 2) || (ny >= mesh.y / 2) || (nz >= mesh.z / 2))
    return;

  int lind = ((nx + mesh.x / 2) * mesh.y * mesh.z + (ny + mesh.y / 2) * mesh.z +
              (nz + mesh.z / 2));

  double alpha_L_i = 1. / alpha_L;
  double n2;
  double U2, ex, ex2;
  int nmx, nmy, nmz;
  double alias1, alias2, alias3, alias4, alias5, alias6;

  alias1 = alias2 = alias3 = alias4 = alias5 = alias6 = 0;

  if ((nx != 0) && (ny != 0) && (nz != 0)) {
    for (int mx = -1; mx <= 1; mx++) {
      nmx = nx + mx * mesh.x;
      for (int my = -1; my <= 1; my++) {
        nmy = ny + my * mesh.y;
        for (int mz = -1; mz <= 1; mz++) {
          nmz = nz + mz * mesh.z;

          n2 = sqr(nmx) + sqr(nmy) + sqr(nmz);

          ex = exp(-(PI * alpha_L_i) * (PI * alpha_L_i) * n2);

          ex2 = sqr(ex);

          U2 = pow((double)csinc(meshi.x * nmx) * csinc(meshi.y * nmy) *
                       csinc(meshi.z * nmz),
                   2.0 * cao);

          alias1 += ex2 / n2;
          alias2 += U2 * ex;
          alias3 += U2 * n2;
          alias4 += U2;

          if (((mx + my + mz) % 2) == 0) { // even term
            alias5 += U2 * n2;
            alias6 += U2;
          } else { // odd term: minus sign!
            alias5 -= U2 * n2;
            alias6 -= U2;
          }
        }
      }
    }

    he_q[lind] =
        (alias1 - sqr(alias2) / (0.5 * (alias3 * alias4 + alias5 * alias6)));

  } else {
    he_q[lind] = 0;
  }
}

double p3m_k_space_error_gpu(double prefactor, int *mesh, int cao, int npart,
                             double sum_q2, double alpha_L, double *box) {
  static thrust::device_vector<double> he_q;

  const size_t mesh_size = mesh[0] * mesh[1] * mesh[2];

  he_q.resize(mesh_size);

  dim3 grid(std::max<int>(1, mesh[0] / 8 + 1),
            std::max<int>(1, mesh[1] / 8 + 1),
            std::max<int>(1, mesh[2] / 8 + 1));

  dim3 block(8, 8, 8);

  int3 mesh3;
  mesh3.x = mesh[0];
  mesh3.y = mesh[1];
  mesh3.z = mesh[2];

  double3 meshi;
  meshi.x = 1. / mesh[0];
  meshi.y = 1. / mesh[1];
  meshi.z = 1. / mesh[2];

  switch (cao) {
  case 1:
    KERNELCALL(p3m_k_space_error_gpu_kernel_ik<1>, grid, block,
               (mesh3, meshi, alpha_L, thrust::raw_pointer_cast(he_q.data())));
    break;
  case 2:
    KERNELCALL(p3m_k_space_error_gpu_kernel_ik<2>, grid, block,
               (mesh3, meshi, alpha_L, thrust::raw_pointer_cast(he_q.data())));
    break;
  case 3:
    KERNELCALL(p3m_k_space_error_gpu_kernel_ik<3>, grid, block,
               (mesh3, meshi, alpha_L, thrust::raw_pointer_cast(he_q.data())));
    break;
  case 4:
    KERNELCALL(p3m_k_space_error_gpu_kernel_ik<4>, grid, block,
               (mesh3, meshi, alpha_L, thrust::raw_pointer_cast(he_q.data())));
    break;
  case 5:
    KERNELCALL(p3m_k_space_error_gpu_kernel_ik<5>, grid, block,
               (mesh3, meshi, alpha_L, thrust::raw_pointer_cast(he_q.data())));
    break;
  case 6:
    KERNELCALL(p3m_k_space_error_gpu_kernel_ik<6>, grid, block,
               (mesh3, meshi, alpha_L, thrust::raw_pointer_cast(he_q.data())));
    break;
  case 7:
    KERNELCALL(p3m_k_space_error_gpu_kernel_ik<7>, grid, block,
               (mesh3, meshi, alpha_L, thrust::raw_pointer_cast(he_q.data())));
    break;
  }

  auto const he_q_final = thrust::reduce(he_q.begin(), he_q.end());

  return 2.0 * prefactor * sum_q2 * sqrt(he_q_final / npart) /
         (box[1] * box[2]);
}
