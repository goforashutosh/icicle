// WARNING: This file is auto-generated by a script.
// Any changes made to this file may be overwritten.
// Please modify the code generation script instead.
// Path to the code generation script: scripts/gen_c_api.py

#pragma once
#ifndef BW6_761_API_H
#define BW6_761_API_H

#include <cuda_runtime.h>
#include "gpu-utils/device_context.cuh"
#include "merkle-tree/merkle.cuh"
#include "curves/params/bw6_761.cuh"
#include "ntt/ntt.cuh"
#include "msm/msm.cuh"
#include "vec_ops/vec_ops.cuh"
#include "poseidon/poseidon.cuh"

extern "C" cudaError_t bw6_761_g2_precompute_msm_bases_cuda(
  bw6_761::g2_affine_t* bases,
  int bases_size,
  int precompute_factor,
  int _c,
  bool are_bases_on_device,
  device_context::DeviceContext& ctx,
  bw6_761::g2_affine_t* output_bases);

extern "C" cudaError_t bw6_761_g2_msm_cuda(
  const bw6_761::scalar_t* scalars, const bw6_761::g2_affine_t* points, int msm_size, msm::MSMConfig& config, bw6_761::g2_projective_t* out);

extern "C" cudaError_t bw6_761_precompute_msm_bases_cuda(
  bw6_761::affine_t* bases,
  int bases_size,
  int precompute_factor,
  int _c,
  bool are_bases_on_device,
  device_context::DeviceContext& ctx,
  bw6_761::affine_t* output_bases);

extern "C" cudaError_t bw6_761_msm_cuda(
  const bw6_761::scalar_t* scalars, const bw6_761::affine_t* points, int msm_size, msm::MSMConfig& config, bw6_761::projective_t* out);

extern "C" bool bw6_761_g2_eq(bw6_761::g2_projective_t* point1, bw6_761::g2_projective_t* point2);

extern "C" void bw6_761_g2_to_affine(bw6_761::g2_projective_t* point, bw6_761::g2_affine_t* point_out);

extern "C" void bw6_761_g2_generate_projective_points(bw6_761::g2_projective_t* points, int size);

extern "C" void bw6_761_g2_generate_affine_points(bw6_761::g2_affine_t* points, int size);

extern "C" cudaError_t bw6_761_g2_affine_convert_montgomery(
  bw6_761::g2_affine_t* d_inout, size_t n, bool is_into, device_context::DeviceContext& ctx);

extern "C" cudaError_t bw6_761_g2_projective_convert_montgomery(
  bw6_761::g2_projective_t* d_inout, size_t n, bool is_into, device_context::DeviceContext& ctx);

extern "C" cudaError_t bw6_761_ecntt_cuda(
  const bw6_761::projective_t* input, int size, ntt::NTTDir dir, ntt::NTTConfig<bw6_761::scalar_t>& config, bw6_761::projective_t* output);

extern "C" bool bw6_761_eq(bw6_761::projective_t* point1, bw6_761::projective_t* point2);

extern "C" void bw6_761_to_affine(bw6_761::projective_t* point, bw6_761::affine_t* point_out);

extern "C" void bw6_761_generate_projective_points(bw6_761::projective_t* points, int size);

extern "C" void bw6_761_generate_affine_points(bw6_761::affine_t* points, int size);

extern "C" cudaError_t bw6_761_affine_convert_montgomery(
  bw6_761::affine_t* d_inout, size_t n, bool is_into, device_context::DeviceContext& ctx);

extern "C" cudaError_t bw6_761_projective_convert_montgomery(
  bw6_761::projective_t* d_inout, size_t n, bool is_into, device_context::DeviceContext& ctx);

extern "C" cudaError_t bw6_761_poseidon_create_cuda(
  poseidon::Poseidon<bw6_761::scalar_t>** poseidon,
  unsigned int arity,
  unsigned int alpha,
  unsigned int partial_rounds,
  unsigned int full_rounds_half,
  const bw6_761::scalar_t* round_constants,
  const bw6_761::scalar_t* mds_matrix,
  const bw6_761::scalar_t* non_sparse_matrix,
  const bw6_761::scalar_t* sparse_matrices,
  const bw6_761::scalar_t domain_tag,
  device_context::DeviceContext& ctx);

extern "C" cudaError_t bw6_761_poseidon_load_cuda(
  poseidon::Poseidon<bw6_761::scalar_t>** poseidon,
  unsigned int arity,
  device_context::DeviceContext& ctx);

extern "C" cudaError_t bw6_761_poseidon_absorb_many_cuda(
  const poseidon::Poseidon<bw6_761::scalar_t>* poseidon,
  const bw6_761::scalar_t* inputs,
  bw6_761::scalar_t* states,
  unsigned int number_of_states,
  unsigned int input_block_len,
  hash::SpongeConfig& cfg);

extern "C" cudaError_t bw6_761_poseidon_squeeze_many_cuda(
  const poseidon::Poseidon<bw6_761::scalar_t>* poseidon,
  const bw6_761::scalar_t* states,
  bw6_761::scalar_t* output,
  unsigned int number_of_states,
  unsigned int output_len,
  hash::SpongeConfig& cfg);

extern "C" cudaError_t bw6_761_poseidon_hash_many_cuda(
  const poseidon::Poseidon<bw6_761::scalar_t>* poseidon,
  const bw6_761::scalar_t* inputs,
  bw6_761::scalar_t* output,
  unsigned int number_of_states,
  unsigned int input_block_len,
  unsigned int output_len,
  hash::SpongeConfig& cfg);

extern "C" cudaError_t
  bw6_761_poseidon_delete_cuda(poseidon::Poseidon<bw6_761::scalar_t>* poseidon, device_context::DeviceContext& ctx);

extern "C" cudaError_t bw6_761_build_poseidon_merkle_tree(
  const bw6_761::scalar_t* leaves,
  bw6_761::scalar_t* digests,
  unsigned int height,
  unsigned int input_block_len, 
  const poseidon::Poseidon<bw6_761::scalar_t>* poseidon_compression,
  const poseidon::Poseidon<bw6_761::scalar_t>* poseidon_sponge,
  const merkle_tree::TreeBuilderConfig& tree_config);

extern "C" cudaError_t bw6_761_mul_cuda(
  bw6_761::scalar_t* vec_a, bw6_761::scalar_t* vec_b, int n, vec_ops::VecOpsConfig& config, bw6_761::scalar_t* result);

extern "C" cudaError_t bw6_761_add_cuda(
  bw6_761::scalar_t* vec_a, bw6_761::scalar_t* vec_b, int n, vec_ops::VecOpsConfig& config, bw6_761::scalar_t* result);

extern "C" cudaError_t bw6_761_sub_cuda(
  bw6_761::scalar_t* vec_a, bw6_761::scalar_t* vec_b, int n, vec_ops::VecOpsConfig& config, bw6_761::scalar_t* result);

extern "C" cudaError_t bw6_761_transpose_matrix_cuda(
  const bw6_761::scalar_t* input,
  uint32_t row_size,
  uint32_t column_size,
  bw6_761::scalar_t* output,
  device_context::DeviceContext& ctx,
  bool on_device,
  bool is_async);

extern "C" cudaError_t bw6_761_bit_reverse_cuda(
  const bw6_761::scalar_t* input,
  uint64_t n,
  vec_ops::BitReverseConfig& config,
  bw6_761::scalar_t* output);

extern "C" void bw6_761_generate_scalars(bw6_761::scalar_t* scalars, int size);

extern "C" cudaError_t bw6_761_scalar_convert_montgomery(
  bw6_761::scalar_t* d_inout, size_t n, bool is_into, device_context::DeviceContext& ctx);

extern "C" cudaError_t bw6_761_initialize_domain(
  bw6_761::scalar_t* primitive_root, device_context::DeviceContext& ctx, bool fast_twiddles_mode);

extern "C" cudaError_t bw6_761_ntt_cuda(
  const bw6_761::scalar_t* input, int size, ntt::NTTDir dir, ntt::NTTConfig<bw6_761::scalar_t>& config, bw6_761::scalar_t* output);

extern "C" cudaError_t bw6_761_release_domain(device_context::DeviceContext& ctx);

#endif