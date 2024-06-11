// WARNING: This file is auto-generated by a script.
// Any changes made to this file may be overwritten.
// Please modify the code generation script instead.
// Path to the code generation script: scripts/gen_c_api.py

#pragma once
#ifndef BN254_API_H
#define BN254_API_H

#include <cuda_runtime.h>
#include "gpu-utils/device_context.cuh"
#include "merkle-tree/merkle.cuh"
#include "curves/params/bn254.cuh"
#include "ntt/ntt.cuh"
#include "msm/msm.cuh"
#include "vec_ops/vec_ops.cuh"
#include "poseidon/poseidon.cuh"
#include "poseidon2/poseidon2.cuh"

extern "C" cudaError_t bn254_g2_precompute_msm_bases_cuda(
  bn254::g2_affine_t* bases,
  int bases_size,
  int precompute_factor,
  int _c,
  bool are_bases_on_device,
  device_context::DeviceContext& ctx,
  bn254::g2_affine_t* output_bases);

extern "C" cudaError_t bn254_g2_msm_cuda(
  const bn254::scalar_t* scalars, const bn254::g2_affine_t* points, int msm_size, msm::MSMConfig& config, bn254::g2_projective_t* out);

extern "C" cudaError_t bn254_precompute_msm_bases_cuda(
  bn254::affine_t* bases,
  int bases_size,
  int precompute_factor,
  int _c,
  bool are_bases_on_device,
  device_context::DeviceContext& ctx,
  bn254::affine_t* output_bases);

extern "C" cudaError_t bn254_msm_cuda(
  const bn254::scalar_t* scalars, const bn254::affine_t* points, int msm_size, msm::MSMConfig& config, bn254::projective_t* out);

extern "C" bool bn254_g2_eq(bn254::g2_projective_t* point1, bn254::g2_projective_t* point2);

extern "C" void bn254_g2_to_affine(bn254::g2_projective_t* point, bn254::g2_affine_t* point_out);

extern "C" void bn254_g2_generate_projective_points(bn254::g2_projective_t* points, int size);

extern "C" void bn254_g2_generate_affine_points(bn254::g2_affine_t* points, int size);

extern "C" cudaError_t bn254_g2_affine_convert_montgomery(
  bn254::g2_affine_t* d_inout, size_t n, bool is_into, device_context::DeviceContext& ctx);

extern "C" cudaError_t bn254_g2_projective_convert_montgomery(
  bn254::g2_projective_t* d_inout, size_t n, bool is_into, device_context::DeviceContext& ctx);

extern "C" cudaError_t bn254_ecntt_cuda(
  const bn254::projective_t* input, int size, ntt::NTTDir dir, ntt::NTTConfig<bn254::scalar_t>& config, bn254::projective_t* output);

extern "C" bool bn254_eq(bn254::projective_t* point1, bn254::projective_t* point2);

extern "C" void bn254_to_affine(bn254::projective_t* point, bn254::affine_t* point_out);

extern "C" void bn254_generate_projective_points(bn254::projective_t* points, int size);

extern "C" void bn254_generate_affine_points(bn254::affine_t* points, int size);

extern "C" cudaError_t bn254_affine_convert_montgomery(
  bn254::affine_t* d_inout, size_t n, bool is_into, device_context::DeviceContext& ctx);

extern "C" cudaError_t bn254_projective_convert_montgomery(
  bn254::projective_t* d_inout, size_t n, bool is_into, device_context::DeviceContext& ctx);

extern "C" cudaError_t bn254_poseidon2_create_cuda(
  poseidon2::Poseidon2<bn254::scalar_t>** poseidon,
  unsigned int width,
  unsigned int alpha,
  unsigned int internal_rounds,
  unsigned int external_rounds,
  const bn254::scalar_t* round_constants,
  const bn254::scalar_t* internal_matrix_diag,
  poseidon2::MdsType mds_type,
  poseidon2::DiffusionStrategy diffusion,
  device_context::DeviceContext& ctx
);

extern "C" cudaError_t bn254_poseidon2_load_cuda(
  poseidon2::Poseidon2<bn254::scalar_t>** poseidon,
  unsigned int width,
  poseidon2::MdsType mds_type,
  poseidon2::DiffusionStrategy diffusion,
  device_context::DeviceContext& ctx
);

extern "C" cudaError_t bn254_poseidon2_absorb_many_cuda(
  const poseidon2::Poseidon2<bn254::scalar_t>* poseidon,
  const bn254::scalar_t* inputs,
  bn254::scalar_t* states,
  unsigned int number_of_states,
  unsigned int input_block_len,
  hash::SpongeConfig& cfg);

extern "C" cudaError_t bn254_poseidon2_squeeze_many_cuda(
  const poseidon2::Poseidon2<bn254::scalar_t>* poseidon,
  const bn254::scalar_t* states,
  bn254::scalar_t* output,
  unsigned int number_of_states,
  unsigned int output_len,
  hash::SpongeConfig& cfg);

extern "C" cudaError_t bn254_poseidon2_hash_many_cuda(
  const poseidon2::Poseidon2<bn254::scalar_t>* poseidon,
  const bn254::scalar_t* inputs,
  bn254::scalar_t* output,
  unsigned int number_of_states,
  unsigned int input_block_len,
  unsigned int output_len,
  hash::SpongeConfig& cfg);

extern "C" cudaError_t
  bn254_poseidon2_delete_cuda(poseidon2::Poseidon2<bn254::scalar_t>* poseidon, device_context::DeviceContext& ctx);

extern "C" cudaError_t bn254_build_merkle_tree(
  const bn254::scalar_t* leaves,
  bn254::scalar_t* digests,
  unsigned int height,
  unsigned int input_block_len, 
  const hash::SpongeHasher<bn254::scalar_t, bn254::scalar_t>* compression,
  const hash::SpongeHasher<bn254::scalar_t, bn254::scalar_t>* bottom_layer,
  const merkle_tree::TreeBuilderConfig& tree_config);

extern "C" cudaError_t bn254_poseidon_create_cuda(
  poseidon::Poseidon<bn254::scalar_t>** poseidon,
  unsigned int arity,
  unsigned int alpha,
  unsigned int partial_rounds,
  unsigned int full_rounds_half,
  const bn254::scalar_t* round_constants,
  const bn254::scalar_t* mds_matrix,
  const bn254::scalar_t* non_sparse_matrix,
  const bn254::scalar_t* sparse_matrices,
  const bn254::scalar_t domain_tag,
  device_context::DeviceContext& ctx);

extern "C" cudaError_t bn254_poseidon_load_cuda(
  poseidon::Poseidon<bn254::scalar_t>** poseidon,
  unsigned int arity,
  device_context::DeviceContext& ctx);

extern "C" cudaError_t bn254_poseidon_absorb_many_cuda(
  const poseidon::Poseidon<bn254::scalar_t>* poseidon,
  const bn254::scalar_t* inputs,
  bn254::scalar_t* states,
  unsigned int number_of_states,
  unsigned int input_block_len,
  hash::SpongeConfig& cfg);

extern "C" cudaError_t bn254_poseidon_squeeze_many_cuda(
  const poseidon::Poseidon<bn254::scalar_t>* poseidon,
  const bn254::scalar_t* states,
  bn254::scalar_t* output,
  unsigned int number_of_states,
  unsigned int output_len,
  hash::SpongeConfig& cfg);

extern "C" cudaError_t bn254_poseidon_hash_many_cuda(
  const poseidon::Poseidon<bn254::scalar_t>* poseidon,
  const bn254::scalar_t* inputs,
  bn254::scalar_t* output,
  unsigned int number_of_states,
  unsigned int input_block_len,
  unsigned int output_len,
  hash::SpongeConfig& cfg);

extern "C" cudaError_t
  bn254_poseidon_delete_cuda(poseidon::Poseidon<bn254::scalar_t>* poseidon);

extern "C" cudaError_t bn254_mul_cuda(
  bn254::scalar_t* vec_a, bn254::scalar_t* vec_b, int n, vec_ops::VecOpsConfig& config, bn254::scalar_t* result);

extern "C" cudaError_t bn254_add_cuda(
  bn254::scalar_t* vec_a, bn254::scalar_t* vec_b, int n, vec_ops::VecOpsConfig& config, bn254::scalar_t* result);

extern "C" cudaError_t bn254_sub_cuda(
  bn254::scalar_t* vec_a, bn254::scalar_t* vec_b, int n, vec_ops::VecOpsConfig& config, bn254::scalar_t* result);

extern "C" cudaError_t bn254_transpose_matrix_cuda(
  const bn254::scalar_t* input,
  uint32_t row_size,
  uint32_t column_size,
  bn254::scalar_t* output,
  device_context::DeviceContext& ctx,
  bool on_device,
  bool is_async);

extern "C" cudaError_t bn254_bit_reverse_cuda(
  const bn254::scalar_t* input,
  uint64_t n,
  vec_ops::BitReverseConfig& config,
  bn254::scalar_t* output);

extern "C" void bn254_generate_scalars(bn254::scalar_t* scalars, int size);

extern "C" cudaError_t bn254_scalar_convert_montgomery(
  bn254::scalar_t* d_inout, size_t n, bool is_into, device_context::DeviceContext& ctx);

extern "C" cudaError_t bn254_initialize_domain(
  bn254::scalar_t* primitive_root, device_context::DeviceContext& ctx, bool fast_twiddles_mode);

extern "C" cudaError_t bn254_ntt_cuda(
  const bn254::scalar_t* input, int size, ntt::NTTDir dir, ntt::NTTConfig<bn254::scalar_t>& config, bn254::scalar_t* output);

extern "C" cudaError_t bn254_release_domain(device_context::DeviceContext& ctx);

#endif