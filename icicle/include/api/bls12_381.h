// WARNING: This file is auto-generated by a script.
// Any changes made to this file may be overwritten.
// Please modify the code generation script instead.
// Path to the code generation script: scripts/gen_c_api.py

#pragma once
#ifndef BLS12_381_API_H
#define BLS12_381_API_H

#include <cuda_runtime.h>
#include "gpu-utils/device_context.cuh"
#include "curves/params/bls12_381.cuh"
#include "ntt/ntt.cuh"
#include "msm/msm.cuh"
#include "vec_ops/vec_ops.cuh"
#include "poseidon/poseidon.cuh"
#include "poseidon/tree/merkle.cuh"

extern "C" cudaError_t bls12_381_g2_precompute_msm_bases_cuda(
  bls12_381::g2_affine_t* bases,
  int bases_size,
  int precompute_factor,
  int _c,
  bool are_bases_on_device,
  device_context::DeviceContext& ctx,
  bls12_381::g2_affine_t* output_bases);

extern "C" cudaError_t bls12_381_g2_msm_cuda(
  const bls12_381::scalar_t* scalars, const bls12_381::g2_affine_t* points, int msm_size, msm::MSMConfig& config, bls12_381::g2_projective_t* out);

extern "C" cudaError_t bls12_381_precompute_msm_bases_cuda(
  bls12_381::affine_t* bases,
  int bases_size,
  int precompute_factor,
  int _c,
  bool are_bases_on_device,
  device_context::DeviceContext& ctx,
  bls12_381::affine_t* output_bases);

extern "C" cudaError_t bls12_381_msm_cuda(
  const bls12_381::scalar_t* scalars, const bls12_381::affine_t* points, int msm_size, msm::MSMConfig& config, bls12_381::projective_t* out);

extern "C" bool bls12_381_g2_eq(bls12_381::g2_projective_t* point1, bls12_381::g2_projective_t* point2);

extern "C" void bls12_381_g2_to_affine(bls12_381::g2_projective_t* point, bls12_381::g2_affine_t* point_out);

extern "C" void bls12_381_g2_generate_projective_points(bls12_381::g2_projective_t* points, int size);

extern "C" void bls12_381_g2_generate_affine_points(bls12_381::g2_affine_t* points, int size);

extern "C" cudaError_t bls12_381_g2_affine_convert_montgomery(
  bls12_381::g2_affine_t* d_inout, size_t n, bool is_into, device_context::DeviceContext& ctx);

extern "C" cudaError_t bls12_381_g2_projective_convert_montgomery(
  bls12_381::g2_projective_t* d_inout, size_t n, bool is_into, device_context::DeviceContext& ctx);

extern "C" cudaError_t bls12_381_ecntt_cuda(
  const bls12_381::projective_t* input, int size, ntt::NTTDir dir, ntt::NTTConfig<bls12_381::scalar_t>& config, bls12_381::projective_t* output);

extern "C" bool bls12_381_eq(bls12_381::projective_t* point1, bls12_381::projective_t* point2);

extern "C" void bls12_381_to_affine(bls12_381::projective_t* point, bls12_381::affine_t* point_out);

extern "C" void bls12_381_generate_projective_points(bls12_381::projective_t* points, int size);

extern "C" void bls12_381_generate_affine_points(bls12_381::affine_t* points, int size);

extern "C" cudaError_t bls12_381_affine_convert_montgomery(
  bls12_381::affine_t* d_inout, size_t n, bool is_into, device_context::DeviceContext& ctx);

extern "C" cudaError_t bls12_381_projective_convert_montgomery(
  bls12_381::projective_t* d_inout, size_t n, bool is_into, device_context::DeviceContext& ctx);

extern "C" cudaError_t bls12_381_build_poseidon_merkle_tree(
  const bls12_381::scalar_t* leaves,
  bls12_381::scalar_t* digests,
  unsigned int height,
  unsigned int arity,
  unsigned int input_block_len, 
  const poseidon::Poseidon<bls12_381::scalar_t>* poseidon_compression,
  const poseidon::Poseidon<bls12_381::scalar_t>* poseidon_sponge,
  const hash::SpongeConfig& sponge_config,
  const merkle_tree::TreeBuilderConfig& tree_config);

extern "C" cudaError_t bls12_381_poseidon_create_cuda(
  poseidon::Poseidon<bls12_381::scalar_t>** poseidon,
  unsigned int arity,
  unsigned int alpha,
  unsigned int partial_rounds,
  unsigned int full_rounds_half,
  const bls12_381::scalar_t* round_constants,
  const bls12_381::scalar_t* mds_matrix,
  const bls12_381::scalar_t* non_sparse_matrix,
  const bls12_381::scalar_t* sparse_matrices,
  const bls12_381::scalar_t domain_tag,
  device_context::DeviceContext& ctx);

extern "C" cudaError_t bls12_381_poseidon_load_cuda(
  poseidon::Poseidon<bls12_381::scalar_t>** poseidon,
  unsigned int arity,
  device_context::DeviceContext& ctx);

extern "C" cudaError_t bls12_381_poseidon_permute_many_cuda(
  const poseidon::Poseidon<bls12_381::scalar_t>* poseidon,
  const bls12_381::scalar_t* states,
  bls12_381::scalar_t* output,
  unsigned int number_of_states,
  device_context::DeviceContext& ctx,
  bool is_async);

extern "C" cudaError_t bls12_381_poseidon_compress_many_cuda(
  const poseidon::Poseidon<bls12_381::scalar_t>* poseidon,
  const bls12_381::scalar_t* states,
  bls12_381::scalar_t* output,
  unsigned int number_of_states,
  unsigned int offset,
  bls12_381::scalar_t* perm_output,
  device_context::DeviceContext& ctx,
  bool is_async);

extern "C" cudaError_t
bls12_381_poseidon_delete_cuda(Poseidon<bls12_381::scalar_t>* poseidon, device_context::DeviceContext& ctx);

extern "C" cudaError_t bls12_381_mul_cuda(
  bls12_381::scalar_t* vec_a, bls12_381::scalar_t* vec_b, int n, vec_ops::VecOpsConfig& config, bls12_381::scalar_t* result);

extern "C" cudaError_t bls12_381_add_cuda(
  bls12_381::scalar_t* vec_a, bls12_381::scalar_t* vec_b, int n, vec_ops::VecOpsConfig& config, bls12_381::scalar_t* result);

extern "C" cudaError_t bls12_381_sub_cuda(
  bls12_381::scalar_t* vec_a, bls12_381::scalar_t* vec_b, int n, vec_ops::VecOpsConfig& config, bls12_381::scalar_t* result);

extern "C" cudaError_t bls12_381_transpose_matrix_cuda(
  const bls12_381::scalar_t* input,
  uint32_t row_size,
  uint32_t column_size,
  bls12_381::scalar_t* output,
  device_context::DeviceContext& ctx,
  bool on_device,
  bool is_async);

extern "C" void bls12_381_generate_scalars(bls12_381::scalar_t* scalars, int size);

extern "C" cudaError_t bls12_381_scalar_convert_montgomery(
  bls12_381::scalar_t* d_inout, size_t n, bool is_into, device_context::DeviceContext& ctx);

extern "C" cudaError_t bls12_381_initialize_domain(
  bls12_381::scalar_t* primitive_root, device_context::DeviceContext& ctx, bool fast_twiddles_mode);

extern "C" cudaError_t bls12_381_ntt_cuda(
  const bls12_381::scalar_t* input, int size, ntt::NTTDir dir, ntt::NTTConfig<bls12_381::scalar_t>& config, bls12_381::scalar_t* output);

extern "C" cudaError_t bls12_381_release_domain(device_context::DeviceContext& ctx);

#endif