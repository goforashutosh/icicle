

#define SHMEM_SIZE 64
#define MAX_SHMEM_LOG_SIZE 6

#include "../ntt/kernel_ntt.cu"
// static inline __device__ uint32_t bit_rev(uint32_t num, uint32_t log_size) { return __brev(num) >> (32 - log_size); }

// template <typename S>
// __global__ void inplace_rbo(S* arr, int size){
// 	int tid = blockIdx.x * blockDim.x + threadIdx.x;
// 	S temp = arr[tid];
// 	arr[tid] = arr[bit_rev(tid)];
// 	arr[bit_rev(tid)] = temp;
// }

template <typename S>
__global__ void mult_and_reduce(S *v, S *v_r, S alpha, int nof_results, int jump_size) {
	// Allocate shared memory
	__shared__ S partial_sum[SHMEM_SIZE];

	// Calculate thread ID
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	// Load elements AND do first add of reduction
	// Vector now 2x as long as number of threads, so scale i
	int i = blockIdx.x * (blockDim.x * 2) + threadIdx.x;

	// Store first partial result instead of just the elements
	// partial_sum[threadIdx.x] = v[i] + v[i + blockDim.x];
	// partial_sum[threadIdx.x] = (S::one() - alpha) * v[2*i] + alpha * v[2*i+1] + (S::one() - alpha) * v[2*(i + blockDim.x)] + alpha * v[2*(i + blockDim.x)+1];
	S e1 = v[i] + alpha * (v[i+jump_size] - v[i]);
	S e2 = v[(i + blockDim.x)] + alpha * (v[(i + blockDim.x)+jump_size] - v[i + blockDim.x]);
	// S e1 = v[2*i] + (v[2*i+1] - v[2*i]);
	// S e2 = v[2*(i + blockDim.x)] + (v[2*(i + blockDim.x)+1] - v[2*(i + blockDim.x)]);
	// partial_sum[threadIdx.x] = v[2*i] + v[2*(i + blockDim.x)] + alpha * (v[2*i+1] - v[2*i] + v[2*(i + blockDim.x)+1] - v[2*(i + blockDim.x)]);
	partial_sum[threadIdx.x] = e1 + e2;
	// __syncthreads();
	v[i] = e1;
	v[i + blockDim.x] = e2;
	// for (int j = 0; j < 2; j++)
	// {
	// 	partial_sum[threadIdx.x] = partial_sum[threadIdx.x] * partial_sum[threadIdx.x];
	// }
	
	__syncthreads();

	// Start at 1/2 block stride and divide by two each iteration
	for (int s = blockDim.x / 2; s > 0; s >>= 1) {
	// for (int s = blockDim.x / 2; s > 1; s >>= 1) {
		// Each thread does work unless it is further than the stride
		if (threadIdx.x < s) {
			partial_sum[threadIdx.x] = partial_sum[threadIdx.x] + partial_sum[threadIdx.x + s];
		}
		__syncthreads();
	}

	// Let the thread 0 for this block write it's result to main memory
	// Result is inexed by this block
	// if (threadIdx.x < nof_results) {
	if (threadIdx.x == 0) {
		// printf("debug tid %d, val %d\n", threadIdx.x, partial_sum[threadIdx.x]);
		// v_r[nof_results*blockIdx.x + threadIdx.x] = partial_sum[threadIdx.x];
		v_r[blockIdx.x] = partial_sum[0];
	}
}


template <typename S>
__global__ void sum_reduction(S *v, S *v_r, int nof_results) {
	// Allocate shared memory
	__shared__ S partial_sum[SHMEM_SIZE];

	// Calculate thread ID
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	// Load elements AND do first add of reduction
	// Vector now 2x as long as number of threads, so scale i
	int i = blockIdx.x * (blockDim.x * 2) + threadIdx.x;

	// Store first partial result instead of just the elements
	partial_sum[threadIdx.x] = v[i] + v[i + blockDim.x];
	__syncthreads();

	// Start at 1/2 block stride and divide by two each iteration
	for (int s = blockDim.x / 2; s > 0; s >>= 1) {
	// for (int s = blockDim.x / 2; s > 1; s >>= 1) {
		// Each thread does work unless it is further than the stride
		if (threadIdx.x < s) {
			partial_sum[threadIdx.x] = partial_sum[threadIdx.x] + partial_sum[threadIdx.x + s];
		}
		__syncthreads();
	}

	// Let the thread 0 for this block write it's result to main memory
	// Result is inexed by this block
	// if (threadIdx.x < nof_results) {
	if (threadIdx.x == 0) {
		// printf("debug tid %d, val %d\n", threadIdx.x, partial_sum[threadIdx.x]);
		v_r[blockIdx.x] = partial_sum[0];
		// v_r[nof_results*blockIdx.x + threadIdx.x] = partial_sum[threadIdx.x];
	}
}

template <typename S>
__global__ void update_evals_kernel(S* evals, S alpha, int poly_size, int nof_ploys){
  int threads_per_poly = poly_size/2;
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= threads_per_poly*nof_ploys) return;
	int poly_id = tid / threads_per_poly;
	int eval_id = tid % threads_per_poly;
  // evals[tid] = (S::one() - alpha) * evals[2*tid] + alpha * evals[2*tid+1];
  // evals[tid] =  evals[2*tid] + (evals[2*tid+1] - evals[2*tid]);
	// if (tid==0) printf("%d, %d, %d, %d, %d\n", poly_size, poly_id, eval_id, poly_id*poly_size*2+eval_id, poly_id*poly_size*2+eval_id+poly_size);
	// if (tid==0) printf("what12 %d %d\n",evals[poly_id*poly_size*2 + eval_id], evals[poly_id*poly_size*2 + eval_id+poly_size]);
  evals[poly_id*poly_size + eval_id] =  evals[poly_id*poly_size+eval_id] + alpha * (evals[poly_id*poly_size+eval_id+threads_per_poly] - evals[poly_id*poly_size+eval_id]);
	// if (tid==0) printf("what %d\n",evals[poly_id*poly_size*2 + eval_id]);
  // evals[tid] = (1 - alpha) * evals[2*tid] + alpha * evals[2*tid+1];
}

template <typename S>
void accumulate(S* in, S* out, int log_size, int nof_results, cudaStream_t stream){
  int nof_steps = (log_size - 1) / MAX_SHMEM_LOG_SIZE;
  int last_step_size = (log_size - 1) % MAX_SHMEM_LOG_SIZE;
	// printf("a nof steps %d last size %d\n", nof_steps, last_step_size);
  for (int i = 0; i < nof_steps; i++)
  {
    sum_reduction<<<1<<(log_size -(MAX_SHMEM_LOG_SIZE)*(i+1)), SHMEM_SIZE/2,0,stream>>>(i? out : in, out, 1);
		// printf("a nof blocks %d\n", 1<<(log_size -(MAX_SHMEM_LOG_SIZE)*(i+1)));
		// cudaDeviceSynchronize();
  	// printf("cuda err %d\n", cudaGetLastError());
  }
  if (last_step_size) sum_reduction<<<2, 1<<(last_step_size-1), 0,stream>>>(nof_steps? out : in, out, 1);
	// cudaDeviceSynchronize();
  // printf("cuda err last %d\n", cudaGetLastError());
}

template <typename S>
void mult_and_accumulate(S* in, S* out, int log_size, S alpha, int nof_results, cudaStream_t stream){
  int nof_steps = (log_size - 1) / MAX_SHMEM_LOG_SIZE;
  int last_step_size = (log_size - 1) % MAX_SHMEM_LOG_SIZE;
	// printf("m nof steps %d last size %d\n", nof_steps, last_step_size);
  for (int i = 0; i < nof_steps; i++)
  {
		if (i) sum_reduction<<<1<<(log_size -(MAX_SHMEM_LOG_SIZE)*(i+1)), SHMEM_SIZE/2,0,stream>>>(i? out : in, out, nof_results);
    else mult_and_reduce<<<1<<(log_size -(MAX_SHMEM_LOG_SIZE)*(i+1)), SHMEM_SIZE/2,0,stream>>>(i? out : in, out, alpha, nof_results, 1<<log_size);
		// if (i) printf("r nof blocks %d\n", 1<<(log_size-(MAX_SHMEM_LOG_SIZE)*(i+1)));
		// else printf("m nof blocks %d\n", 1<<(log_size-(MAX_SHMEM_LOG_SIZE)*(i+1)));
		// cudaDeviceSynchronize();
  	// printf("cuda err %d\n", cudaGetLastError());
  }
  if (last_step_size) {
		if (nof_steps) sum_reduction<<<2, 1<<(last_step_size-1), 0,stream>>>(nof_steps? out : in, out, nof_results);
		else mult_and_reduce<<<2, 1<<(last_step_size-1), 0,stream>>>(nof_steps? out : in, out, alpha, nof_results, 1<<(last_step_size+1));
		// if (nof_steps) printf("r last");
		// else printf("m last");
	} 
	// cudaDeviceSynchronize();
  // printf("cuda err last %d\n", cudaGetLastError());
}

template <typename S>
 __launch_bounds__(1)
__global__ void add_to_trace(S* trace, S* vals, int n, int round_num, int nof_polys){
	// for (int i = 0; i < nof_polys+1; i++)
	// {
	// 	trace[2*round_num+1+i] = vals[i];
	// }
	  trace[2*round_num+1] = vals[0];
    trace[2*round_num+2] = vals[1];
		// printf("%d  %d\n", vals[0], vals[1]);
}

template <typename S>
// __global__ void combinations_kernel(S* in, S* out, S (*combine_func)()){
__global__ void combinations_kernel3(S* in, S* out, int poly_size){
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	S rp[4] = {S::one(), S::one(), S::one(), S::one()};
	S e1, e2;
	#pragma unroll
	for (int l = 0; l < 3; l++)
	{
	  e1 = in[l*poly_size + 2*tid];
	  e2 = in[l*poly_size + 2*tid + 1];
		rp[0] = l? rp[0]*e1 : e1; //k=0
		rp[1] = l? rp[1]*e2 : e2; //k=1
		rp[2] = l? rp[2]*(e2 + e2 - e1) : (e2 + e2 - e1); //k=2
		rp[3] = l? rp[3]*(e1 + e1 - e2) : (e1 + e1 - e2); //k=-1
		// rp[3] = l? rp[3]*(e2 + e2 + e2 - e1 - e1) : (e2 + e2 + e2 - e1 - e1); //k=3
	}
	out[4*tid] = rp[0];
	out[4*tid+1] = rp[1];
	out[4*tid+2] = rp[2];
	out[4*tid+3] = rp[3];
}

template <typename S>
// __global__ void combinations_kernel(S* in, S* out, S (*combine_func)()){
__global__ void mult_and_combine3(S* in, S* out, int poly_size, S alpha){
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	S rp[4] = {S::one(), S::one(), S::one(), S::one()};
	S e1, e2;
	#pragma unroll
	for (int l = 0; l < 3; l++)
	{
		e1 = in[l*poly_size + 4*tid] + alpha * (in[l*poly_size + 4*tid + 1] - in[l*poly_size + 4*tid]);
		e2 = in[l*poly_size + 4*tid + 2] + alpha * (in[l*poly_size + 4*tid + 3] - in[l*poly_size + 4*tid + 2]);
		rp[0] = rp[0]*e1;
		rp[1] = rp[1]*e2;
		rp[2] = rp[2]*(e2 + e2 - e1);
		rp[3] = rp[3]*(e2 + e2 + e2 - e1 - e1);
	}
	out[4*tid] = rp[0];
	out[4*tid+1] = rp[1];
	out[4*tid+2] = rp[2];
	out[4*tid+3] = rp[3];
}

// template <typename S, int M>
// // __global__ void combinations_kernel(S* in, S* out, S (*combine_func)()){
// __global__ void combinations_kernel(S* in, S* out){
// 	int tid = blockIdx.x * blockDim.x + threadIdx.x;
// 	S rp = s::one;
// 	#pragma unroll
// 	for (int k = 0; k < M+1; k++)
// 	{
// 		#pragma unroll
// 		for (int l = 0; l < M; i++)
// 		{
// 			rp *= in[2*tid] * (1 - k) + in[2*tid + 1] * k;
// 		}
// 	}
// }

// template <typename S>
// __device__ S simple_combine(S* f){
// 	return f[0]*f[1]*f[2]
// }

template <typename S>
S my_hash(){
	S val = S::one() + S::one();
	val = val + val;
	val = val + val; 
	return val + S::one() + S::one();
}

template <typename S>
void sumcheck_alg1(S* evals, S* t, S* T, S C, int n, bool reorder, cudaStream_t stream){
	if (reorder) reorder_digits_inplace_and_normalize_kernel<<<1<<(max(n-6,0)),64,0,stream>>>(evals, n, false, ntt::eRevType::NaturalToRev, false, S::one());
	// S alpha = 1;
	// S alpha = S::one();
	S alpha = my_hash<S>();
	// S alpha = S::rand_host();
  // S alpha = my_hash(/*T, C*/);
  // S rp_even, rp_odd;
  for (int p = 0; p < n-1; p++)
  {
    int nof_threads = 1<<(n-1-p);
		printf("reg nof threads %d\n", nof_threads);
    // move update kernel here and unify
    // reduction_kernel<<<nof_threads>>>(evals, t, n-p); //accumulation
    accumulate(evals, t, n-p, 2, stream); //accumulation
		// cudaDeviceSynchronize();
		// printf("cuda a err %d\n", cudaGetLastError());
		add_to_trace<<<1,1,0,stream>>>(T, t, n, p, 1);
		// cudaDeviceSynchronize();
		// printf("cuda t err %d\n", cudaGetLastError());
    // T[2*p+1] = t[0];
    // T[2*p+2] = t[1];
    // alpha = my_hash(/*alpha, t[0], t[1]*/); //phase 2
		int NOF_THREADS = min(256,nof_threads);
		int NOF_BLOCKS = (nof_threads + NOF_THREADS - 1) / NOF_THREADS;
    update_evals_kernel<<<NOF_BLOCKS, NOF_THREADS,0, stream>>>(evals, alpha, 1<<(n-p), 1); //phase 3
		// cudaDeviceSynchronize();
		// printf("cuda err u %d\n", cudaGetLastError());
		#ifdef DEBUG
		break;
		#endif
  }
	add_to_trace<<<1,1,0,stream>>>(T, evals, n, n-1, 1);
}

template <typename S>
void sumcheck_alg1_unified(S* evals, S* t, S* T, S C, int n, bool reorder, cudaStream_t stream){
	if (reorder) reorder_digits_inplace_and_normalize_kernel<<<1<<(max(n-6,0)),64,0,stream>>>(evals, n, false, ntt::eRevType::NaturalToRev, false, S::one());
	// S alpha = 1;
	// S alpha = S::one() + S::one();
	S alpha = my_hash<S>();
	// S alpha = S::rand_host();
  // S alpha = my_hash(/*T, C*/);
  // S rp_even, rp_odd;
  for (int p = 0; p < n-1; p++)
  // for (int p = 0; p < 2; p++)
  {
    int nof_threads = 1<<(n-1-p);
		// printf("nof threads %d\n", nof_threads);
    // move update kernel here and unify
    // reduction_kernel<<<nof_threads>>>(evals, t, n-p); //accumulation
    if (p) mult_and_accumulate(evals, t, n-p, alpha, 2, stream); //accumulation
		else accumulate(evals, t, n-p, 2, stream);
		add_to_trace<<<1,1,0,stream>>>(T, t, n, p, 1);
    // T[2*p+1] = t[0];
    // T[2*p+2] = t[1];
    // alpha = my_hash(/*alpha, t[0], t[1]*/); //phase 2
		// int NOF_THREADS = 256;
		// int NOF_BLOCKS = (nof_threads + NOF_THREADS - 1) / NOF_THREADS;
    // update_evals_kernel<<<NOF_BLOCKS, NOF_THREADS,0, stream>>>(evals, alpha); //phase 3
		#ifdef DEBUG
		if (p) break;
		#endif
  }
	#ifndef DEBUG
	update_evals_kernel<<<1, 2,0, stream>>>(evals, alpha, 4, 1);
	#endif
	add_to_trace<<<1,1,0,stream>>>(T, evals, n, n-1, 1);
}

template <typename S>
void sumcheck_alg3_poly3(S* evals, S* t, S* T, S C, int n, bool reorder, cudaStream_t stream){
	if (reorder) reorder_digits_inplace_and_normalize_kernel<<<1<<(max(n-6,0)),64,0,stream>>>(evals, n, false, ntt::eRevType::NaturalToRev, false, S::one());
	// S alpha = 1;
	// S alpha = S::one();
	S alpha = my_hash<S>();
	// S alpha = S::rand_host();
  // S alpha = my_hash(/*T, C*/);
  // S rp_even, rp_odd;
  for (int p = 0; p < n-1; p++)
  {
    int nof_threads = 1<<(n-1-p);
		int NOF_THREADS = 64;
		int NOF_BLOCKS = (nof_threads + NOF_THREADS - 1) / NOF_THREADS;
		// printf("nof threads %d\n", nof_threads);
    // move update kernel here and unify
    // reduction_kernel<<<nof_threads>>>(evals, t, n-p); //accumulation
		combinations_kernel3<<<NOF_BLOCKS, NOF_THREADS,0,stream>>>(evals, t, 1<<n);
		accumulate(t, t, n-p, 4, stream);
		add_to_trace<<<1,1,0,stream>>>(T, t, p, 3);
    // T[2*p+1] = t[0];
    // T[2*p+2] = t[1];
    // alpha = my_hash(/*alpha, t[0], t[1]*/); //phase 2
    update_evals_kernel<<<NOF_BLOCKS, NOF_THREADS,0, stream>>>(evals, alpha, nof_threads); //phase 3
  }
	// update_evals_kernel<<<1, 2,0, stream>>>(evals, alpha);
	add_to_trace<<<1,1,0,stream>>>(T, evals, n-1, 3);
}

template <typename S>
void sumcheck_alg3_poly3_unified(S* evals, S* t, S* T, S C, int n, cudaStream_t stream){
	// S alpha = 1;
	// S alpha = S::one();
	// S alpha = S::rand_host();
  S alpha = my_hash<S>();
  // S rp_even, rp_odd;
  for (int p = 0; p < n-1; p++)
  {
    int nof_threads = 1<<(n-1-p);
		int NOF_THREADS = 64;
		int NOF_BLOCKS = (nof_threads + NOF_THREADS - 1) / NOF_THREADS;
		// printf("nof threads %d\n", nof_threads);
    // move update kernel here and unify
    // reduction_kernel<<<nof_threads>>>(evals, t, n-p); //accumulation
		if (p) mult_and_combine3<<<NOF_BLOCKS, NOF_THREADS,0,stream>>>(evals, t, 1<<n, alpha);
		else combinations_kernel3<<<NOF_BLOCKS, NOF_THREADS,0,stream>>>(evals, t, 1<<n);
		accumulate(t, t, n-p, 4, stream);
		add_to_trace<<<1,1,0,stream>>>(T, t, p, 3);
    // T[2*p+1] = t[0];
    // T[2*p+2] = t[1];
    // alpha = my_hash(/*alpha, t[0], t[1]*/); //phase 2
    // update_evals_kernel<<<NOF_BLOCKS, NOF_THREADS,0, stream>>>(evals, alpha, nof_threads); //phase 3
  }
	update_evals_kernel<<<1, 2,0, stream>>>(evals, alpha, 2);
	add_to_trace<<<1,1,0,stream>>>(T, evals, n-1, 3);
}


template <typename S>
void sumcheck_alg1_ref(S* evals, S* t, S* T, S C, int n){
  // S alpha = my_hash(/*T, C*/);
	// S alpha = 1;
	// S alpha = S::one() + S::one();
	S alpha = my_hash<S>();
  S rp_bottom, rp_top;
  for (int p = 0; p < n; p++)
  {
		// rp_even = 0; rp_odd = 0;
		rp_bottom = S::zero(); rp_top = S::zero();
		// printf("evals\n");
		// for (int i = 0; i < 1<<(n-p); i++)
		// {
		// 	printf("%d, ",evals[i]);
		// }
		// printf("\n");
		for (int i = 0; i < 1<<(n-1-p); i++)
		{
			rp_bottom = rp_bottom + evals[i];
			rp_top = rp_top + evals[i+(1<<(n-1-p))];
		}
    T[2*p+1] = rp_bottom;
    T[2*p+2] = rp_top;
    // alpha = my_hash(/*alpha, t[0], t[1]*/); //phase 2
		// alpha = 1;
		// alpha = S::one();
		for (int i = 0; i < 1<<(n-1-p); i++)
		{
			t[i] = (S::one() - alpha) * evals[i] + alpha * evals[i+(1<<(n-1-p))];
			// t[i] = (1-alpha)*evals[2*i] + alpha*evals[2*i+1];
		}
		for (int i = 0; i < 1<<(n-1-p); i++)
		{
			evals[i] = t[i];
		}
  }
}