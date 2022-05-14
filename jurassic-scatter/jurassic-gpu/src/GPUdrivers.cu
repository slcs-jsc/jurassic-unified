#include <cuda.h>
#include <omp.h>
#include "jr_common.h" // ...

#ifdef GPUDEBUG
    #define debug_printf(...) printf(__VA_ARGS__)
#else
    #define debug_printf(...)
#endif

	// Helper /////////////////////////////////////////////////////////////////////
	// Checking return types of all CUDA runtime functions is best practice, 
	//  ... has negligible performance impact and should not be omitted unless absolutely necessary
	__host__ inline
	void __cudaSafeCall(cudaError err, const char *file, const int line, char const *call=nullptr) { // Actual check function
		if (cudaSuccess != err) {
			fprintf(stderr, "[ERROR] CUDA call%s%s at %s:%d\n%s\n", call?" to ":"", call, file, line, cudaGetErrorString(err));
			exit(0);
		}
	} // __cudaSafeCall
    #define cuCheck(err) __cudaSafeCall((err), __FILE__, __LINE__, #err) // Syntactic sugar to enhance output

    // As CUDA kernel launches are asynchronous error checking is more difficult, 
    // ... as the check might occur prior to the actual error - this macro makes 
    // ... sure it catches an error if it occurs by explicit Synchronization. 
    // ... Due to the performance impact it is only active in debug mode.
    __host__ inline
    void __cuKernelCheck(const char* file, const int line) {
#ifdef GPUDEBUG
		cudaDeviceSynchronize();
		cudaError_t err = cudaPeekAtLastError();
		if (cudaSuccess != err) {
			fprintf(stderr, "[ERROR] CUDA kernel call at %s:%d\n%s\n",  file, line, cudaGetErrorString(err));
			exit(0);
		} // err
#endif
	} // __cuKernelCheck
    #define cuKernelCheck() __cuKernelCheck(__FILE__, __LINE__)

	// GPU Memory management /////////////////////////////////////////////////////////

    __host__
    void copy_data_to_GPU(void *d, void const *h, size_t const nBytes, cudaStream_t const stream) {
        debug_printf("[INFO] transfer %lu Byte from %p @host to %p @device\n", nBytes, h, d);
        cuCheck(cudaMemcpyAsync(d, h, nBytes, cudaMemcpyHostToDevice, stream));
    } // copy_data_to_GPU

    __host__
    void get_data_from_GPU(void *h, void const *d, size_t const nBytes, cudaStream_t const stream) {
        debug_printf("[INFO] transfer %lu Byte from %p @device to %p @host\n", nBytes, d, h);
        cuCheck(cudaMemcpyAsync(h, d, nBytes, cudaMemcpyDeviceToHost, stream));
    } // get_data_from_GPU

    __host__
    void* __allocate_on_GPU(size_t const nBytes, char const *srcfile=nullptr, int const srcline=0) {
        debug_printf("[INFO] cudaMalloc %.6f MByte in %s:%i\n", 1e-6*nBytes, srcfile, srcline);
        void* d = nullptr;
        cuCheck(cudaMalloc(&d, nBytes));
        return d;
    } // allocate_on_GPU
    #define malloc_GPU(TYPE, NUM) (TYPE *)__allocate_on_GPU((NUM)*sizeof(TYPE), __FILE__, __LINE__)

    __host__
    void free_memory_on_GPU(void**d) {
        cuCheck(cudaFree(*d));
        *d = nullptr;
    } // free_memory_on_GPU

    
    __host__
    void* __allocate_unified_memory(size_t const nBytes, char const *srcfile=nullptr, int const srcline=0) {
        debug_printf("[INFO] cudaMallocManaged %.6f MByte in %s:%i\n", 1e-6*nBytes, srcfile, srcline);
        void* d = nullptr;
        cuCheck(cudaMallocManaged(&d, nBytes));
        return d;
    } // allocate_on_GPU
    #define getUnifiedMemory(TYPE, NUM) (TYPE *)__allocate_unified_memory((NUM)*sizeof(TYPE), __FILE__, __LINE__)
    
  extern "C" {
	  trans_table_t* get_tbl_on_GPU(ctl_t const *ctl); 
  }
    __host__
	trans_table_t* get_tbl_on_GPU(ctl_t const *ctl) {
    printf("DEBUG #%d get_tbl_on_GPU was called..\n", ctl->MPIglobrank);
		static trans_table_t *tbl_G = nullptr;
		if (!tbl_G) {
      printf("DEBUG #%d tbl_G == nullptr\n", ctl->MPIglobrank);
			trans_table_t* tbl = get_tbl(ctl);
#ifdef  USE_UNIFIED_MEMORY_FOR_TABLES
            printf("[INFO] allocated %.3f MByte unified memory for tables\n", 1e-6*sizeof(trans_table_t));
            tbl_G = tbl; // just passing a pointer, same memory space
#else
            printf("[INFO] try to allocate %.3f MByte GPU memory for tables\n", 1e-6*sizeof(trans_table_t));
			tbl_G = malloc_GPU(trans_table_t, 1);
			copy_data_to_GPU(tbl_G, tbl, sizeof(trans_table_t), 0);
#endif
		} // !tbl_G
		return tbl_G;
	} // get_tbl_on_GPU

	// ################ GPU driver routines - keep consistent with CPUdrivers.cu ##############

	// Radiance -> Brightness conversion //////////////////////////////////////////
	void __global__ // GPU-kernel
		radiance_to_brightness_GPU(ctl_t const *ctl, obs_t *obs) { // operates onto obs in-place
			for(int ir = blockIdx.x; ir < obs->nr; ir += gridDim.x) { // grid stride loop over blocks = rays
				for(int id = threadIdx.x; id < ctl->nd; id += blockDim.x) { // grid stride loop over threads = detectors
                    auto const radiance = obs->rad[ir][id];
					obs->rad[ir][id] = brightness_core(radiance, ctl->nu[id]); // modify in-place
				} // id
			} // ir
		} // radiance_to_brightness_GPU

	// Add planetary surface emission ////////////////////////////////////////////
	void __global__ // GPU-kernel
		surface_terms_GPU(const trans_table_t *tbl, obs_t *obs, double const tsurf[], int const nd) {
			for(int ir = blockIdx.x; ir < obs->nr; ir += gridDim.x) { // grid stride loop over blocks = rays
				for(int id = threadIdx.x; id < nd; id += blockDim.x) { // grid stride loop over threads = detectors
					add_surface_core(obs, tbl, tsurf[ir], ir, id);
				} // id
			} // ir
		} // surface_terms_GPU

  template<int CO2, int H2O, int N2, int O2>
  __host__ __device__ __ext_inline__
  double continua_core_temp(ctl_t const *ctl, pos_t const *los, int const ig_co2, int const ig_h2o, int const id) {
      double const p = los->p;
      double const t = los->t;
      double const ds = los->ds;
      double beta_ds = los->k[ctl->window[id]]*ds;													// extinction
      // make sure that ig_co2 and ig_h2o are both >= 0
      if ((CO2)) beta_ds += continua_ctmco2(ctl->nu[id], p, t, los->u[ig_co2]);						// co2 continuum
      if ((H2O)) beta_ds += continua_ctmh2o(ctl->nu[id], p, t, los->q[ig_h2o], los->u[ig_h2o]);		// h2o continuum
      if ((N2))  beta_ds += continua_ctmn2(ctl->nu[id], p, t)*ds;									// n2 continuum
      if ((O2))  beta_ds += continua_ctmo2(ctl->nu[id], p, t)*ds;									// o2 continuum
      return     beta_ds;
  } // continua_core_bbbb where each b is either 0 or 1
  
  template<int CO2, int H2O, int N2, int O2>
    void __global__ // GPU-kernel
    fusion_kernel_GPU(trans_table_t const *tbl, ctl_t const *ctl,
        obs_t *obs, pos_t const (*restrict los)[NLOS],
        int const np[], int const ig_co2, int const ig_h2o,
        double const (*restrict aero_beta)[ND]) {
			double tau_path[NG];
			for(int ir = blockIdx.x; ir < obs->nr; ir += gridDim.x) { // grid stride loop over blocks = rays
				for(int id = threadIdx.x; id < ND; id += blockDim.x) { // grid stride loop over threads = detectors
					obs->rad[ir][id] = 0.0;
					obs->tau[ir][id] = 1.0;
				} // id
				for(int ig = 0; ig < NG; ++ig) {
					tau_path[ig] = 1.0;
				} // ig
				for(int ip = 0; ip < np[ir]; ++ip) {
					for(int id = threadIdx.x; id < ctl->nd; id += blockDim.x) { // grid stride loop over threads = detectors
						double const beta_ds = continua_core_temp <CO2, H2O, N2, O2>
                                   (ctl, &(los[ir][ip]), ig_co2, ig_h2o, id);           // function args
						
            //Added:
            double aero_ds = 0;
            if(NULL != aero_beta) // only if scattering is included
              aero_ds = los[ir][ip].aerofac * aero_beta[los[ir][ip].aeroi][id] * los[ir][ip].ds;
            
            double const tau_gas = apply_ega_core(tbl, &(los[ir][ip]), tau_path, ctl->ng, id);
						double const planck = src_planck_core(tbl, los[ir][ip].t, id);
						new_obs_core(obs, ir, id, beta_ds + aero_ds, planck, tau_gas);
					} // id --> parallel over detectors=threads
				} // ip --> non-parallelisable
			} // ir --> parallel over rays==blocks
    } // fusion_kernel_GPU

    __host__
	void multi_version_GPU(char const fourbit, trans_table_t const *tbl, ctl_t const *ctl,
			obs_t *obs, pos_t const (*restrict los)[NLOS],
			int const np[], int const ig_co2, int const ig_h2o,
      double const (*restrict aero_beta)[ND],
			unsigned const grid, unsigned const block, unsigned const shmem, cudaStream_t const stream) {
#define LaunchKernel <<< grid, block, shmem, stream >>> (tbl, ctl, obs, los, np, ig_co2, ig_h2o, aero_beta)
		switch (fourbit) {
			case 0b0000: fusion_kernel_GPU<0,0,0,0> LaunchKernel; break;
			case 0b0001: fusion_kernel_GPU<0,0,0,1> LaunchKernel; break;
			case 0b0010: fusion_kernel_GPU<0,0,1,0> LaunchKernel; break;
			case 0b0011: fusion_kernel_GPU<0,0,1,1> LaunchKernel; break;
			case 0b0100: fusion_kernel_GPU<0,1,0,0> LaunchKernel; break;
			case 0b0101: fusion_kernel_GPU<0,1,0,1> LaunchKernel; break;
			case 0b0110: fusion_kernel_GPU<0,1,1,0> LaunchKernel; break;
			case 0b0111: fusion_kernel_GPU<0,1,1,1> LaunchKernel; break;
			case 0b1000: fusion_kernel_GPU<1,0,0,0> LaunchKernel; break;
			case 0b1001: fusion_kernel_GPU<1,0,0,1> LaunchKernel; break;
			case 0b1010: fusion_kernel_GPU<1,0,1,0> LaunchKernel; break;
			case 0b1011: fusion_kernel_GPU<1,0,1,1> LaunchKernel; break;
			case 0b1100: fusion_kernel_GPU<1,1,0,0> LaunchKernel; break;
			case 0b1101: fusion_kernel_GPU<1,1,0,1> LaunchKernel; break;
			case 0b1110: fusion_kernel_GPU<1,1,1,0> LaunchKernel; break;
			case 0b1111: fusion_kernel_GPU<1,1,1,1> LaunchKernel; break;
		} // fourbit
#undef	LaunchKernel
	} // multi_version_GPU

	// Raytracing ////////////////////////////////////////////////////////////////
  template<int scattering_included>
	void __global__ // GPU-kernel
		raytrace_rays_GPU(ctl_t const *ctl, const atm_t *atm, obs_t *obs, pos_t
    los[][NLOS], double *tsurf, int np[], int const *atm_id, aero_t const *aero) {
			for(int ir = blockIdx.x*blockDim.x + threadIdx.x; ir < obs->nr; ir += blockDim.x*gridDim.x) { // grid stride loop over rays
        np[ir] = traceray<scattering_included>(ctl, &atm[(NULL == atm_id ? 0 : atm_id[ir])], obs, ir, los[ir], &tsurf[ir], aero);
			} // ir
		} // raytrace_rays_GPU

	// Compute hydrostatic equilibria for all atm //////////////////////////////
	void __global__ // GPU-kernel
		hydrostatic_kernel_GPU(ctl_t const *ctl, atm_t *atm, const int num_of_atms, int const ig_h2o) {
			for(int i = blockIdx.x*blockDim.x + threadIdx.x; i < num_of_atms; i += blockDim.x*gridDim.x) {
				hydrostatic_1d_h2o(ctl, &atm[i], 0, atm[i].np, ig_h2o);
			} // ip
		} // hydrostatic_kernel

    __host__
	void hydrostatic1d_GPU(ctl_t const *ctl, ctl_t const *ctl_G,
			atm_t *atm_G, int const num_of_atms, int const ig_h2o, cudaStream_t const stream) {
		if(ctl->hydz < 0) return; // Check reference height
		hydrostatic_kernel_GPU<<<num_of_atms/32 + 1, 32, 0, stream>>> (ctl_G, atm_G, num_of_atms, ig_h2o);
	} // hydrostatic1d_GPU

	// ################ end of GPU driver routines ##############

	// GPU control struct containing GPU version of input, intermediate and output arrays
	typedef struct {
    aero_t *aero_G; // new
		obs_t  *obs_G;
		atm_t  *atm_G;
		pos_t (*los_G)[NLOS];
		double *tsurf_G;
		int    *np_G;
    double (*aero_beta_G)[ND];
    int *atm_id_G; 
		cudaStream_t stream;
	} gpuLane_t;

	// The full forward model working on one package of NR rays
	void formod_one_package_GPU(ctl_t *ctl, ctl_t *ctl_G,
			trans_table_t const *tbl_G,
			atm_t *atm, // can be made const if we do not get the atms back
			obs_t *obs,
      aero_t const *aero,
      int const *atm_id,
			gpuLane_t const *gpu)
    // a workload manager for the GPU
    {
		debug_printf("[INFO] %s GPU\n"
               " Rays:    %9d (max %d)\n"
               " Gases:   %9d (max %d)\n"
               " Channels:%9d (max %d)\n",
               __func__, obs->nr, NR, ctl->ng, NG, ctl->nd, ND);
        
		atm_t *atm_G = gpu->atm_G;
    int *atm_id_G = gpu->atm_id_G;
    aero_t *aero_G = gpu->aero_G;
		obs_t *obs_G = gpu->obs_G;
    double (*aero_beta_G)[ND] = gpu->aero_beta_G;
    pos_t (* los_G)[NLOS] = gpu->los_G;
		double *tsurf_G = gpu->tsurf_G;
		int *np_G = gpu->np_G;
	  cudaSetDevice(0);
    cudaEvent_t finishedEvent;
		cudaEventCreate(&finishedEvent);

		// gas absorption continua configuration
		static int ig_co2 = -999, ig_h2o = -999;
		if((ctl->ctm_h2o) && (ig_h2o == -999)) ig_h2o = jur_find_emitter(ctl, "H2O");
		if((ctl->ctm_co2) && (ig_co2 == -999)) ig_co2 = jur_find_emitter(ctl, "CO2");
		// binary switches for the four gases
		char const fourbit = (char)
                ( ( (1 == ctl->ctm_co2) && (ig_co2 >= 0) )*0b1000   // CO2
                + ( (1 == ctl->ctm_h2o) && (ig_h2o >= 0) )*0b0100   // H2O
                +   (1 == ctl->ctm_n2)                    *0b0010   // N2
                +   (1 == ctl->ctm_o2)                    *0b0001); // O2

		unsigned const nd = ctl->nd, nr = obs->nr; // abbreviate
    int const num_of_atms = jur_get_num_of_atms(nr, atm_id);

    // "beta_a" -> 'a', "beta_e" -> 'e'
    char const beta_type = ctl->sca_ext[5];
		
    cudaStream_t stream = gpu->stream;

		copy_data_to_GPU(obs_G, obs, 1*sizeof(obs_t), stream);
    if(NULL != aero) { // only if scattering is included
      copy_data_to_GPU(aero_G, aero, 1*sizeof(aero_t), stream);
      copy_data_to_GPU(aero_beta_G, ('a' == beta_type) ? aero->beta_a :
                       aero->beta_e, NLMAX * ND * sizeof(double), stream);
    }
    if(NULL != atm_id) {
      copy_data_to_GPU(atm_id_G, atm_id, nr * sizeof(int), stream);
		  copy_data_to_GPU(atm_G, atm, num_of_atms * sizeof(atm_t), stream);
    }
    else {
      atm_id_G = NULL;
		  copy_data_to_GPU(atm_G, atm, sizeof(atm_t), stream);
    }
        
		hydrostatic1d_GPU(ctl, ctl_G, atm_G, num_of_atms, ig_h2o, stream); // in this call atm_G gets modified
		cuKernelCheck();

    if(NULL != aero && ctl->sca_n > 0) { // only if scattering is included
      raytrace_rays_GPU<1> <<< (nr/64)+1, 64, 0, stream>>> (ctl_G, atm_G, obs_G, los_G, tsurf_G, np_G, atm_id_G, aero_G);
      cuKernelCheck();
    } else {
      raytrace_rays_GPU<0> <<< (nr/64)+1, 64, 0, stream>>> (ctl_G, atm_G, obs_G, los_G, tsurf_G, np_G, atm_id_G, NULL);
      cuKernelCheck();
    }
	
    multi_version_GPU(fourbit, tbl_G, ctl_G, obs_G, los_G, np_G, ig_co2, ig_h2o,
                      NULL == aero ? NULL:  aero_beta_G,
                      nr, nd, ctl->gpu_nbytes_shared_memory, stream);
		cuKernelCheck();

    surface_terms_GPU <<< nr, nd, 0, stream>>> (tbl_G, obs_G, tsurf_G, nd);
		cuKernelCheck();
        
    if (ctl->write_bbt && ctl->leaf_nr == -1) { // convert radiance to brightness (in-place)
        radiance_to_brightness_GPU <<< nr, nd, 0, stream >>> (ctl_G, obs_G);
    } // write_bbt

// 		get_data_from_GPU(atm, atm_G, 1*sizeof(atm_t), stream); // do we really need to get the atms back?
		get_data_from_GPU(obs, obs_G, 1*sizeof(obs_t), stream); // always transfer NR rays
    
    printf("Heisenbug :) %lf\n", obs->rad[0][0]);

		// Wait for GPU operations to complete
		cuCheck(cudaEventRecord(finishedEvent, stream));
		cuCheck(cudaEventSynchronize(finishedEvent));
	} // formod_one_package_GPU

    // make sure that jur_formod_multiple_packages_GPU can be linked from CPUdrivers.c
	extern "C" {
	   void jur_formod_multiple_packages_GPU(ctl_t *ctl, atm_t *atm, obs_t *obs,
                                           int n, int const *atm_id, aero_t const *aero);
   }

	__host__
	void jur_formod_multiple_packages_GPU(ctl_t *ctl, atm_t *atm, obs_t *obs,
                                        int n, int const *atm_id, aero_t const *aero) {
    static ctl_t *ctl_G=NULL;
		static trans_table_t *tbl_G=NULL;

		static int numDevices = 0;
		static gpuLane_t* gpuLanes=NULL;
		static size_t numLanes = 0;
		// static size_t nextLane = 0;
		// size_t myLane = 0;

		static bool do_init = true;

    // it can also by set to true, the initial value is irelevant
    static bool multi_atm_before = false;
    const bool multi_atm_now = NULL != atm_id;

    atm_t **divided_atms;
    int **divided_atm_ids;

    if(multi_atm_now) {
      divided_atms = (atm_t **) malloc(n * sizeof(atm_t *));
      divided_atm_ids = (int **) malloc(n * sizeof(int *));
      jur_divide_atm_data_into_packages(atm, obs, n, atm_id, divided_atms, divided_atm_ids);
    }

#pragma omp critical
		{
			if (do_init || (!do_init && multi_atm_before != multi_atm_now)) {
              double tic = omp_get_wtime();
              size_t sizePerLane = sizeof(aero_t) + sizeof(obs_t) + sizeof(atm_t) + NR * (NLOS * sizeof(pos_t) + sizeof(double) + sizeof(int));
              // in this case we have NR * atm_t instead of the only one and one additional atm_id array
              if(multi_atm_now)
                sizePerLane = sizeof(aero_t) + sizeof(obs_t) + NR * (sizeof(atm_t) + NLOS * sizeof(pos_t) + sizeof(double) + 2 * sizeof(int));
 
              if (ctl->checkmode) {
                printf("# %s: GPU memory requirement per lane is %.3f MByte\n", __func__, 1e-6*sizePerLane);
              } else {
                cudaSetDevice(0);
                cuCheck(cudaGetDeviceCount(&numDevices));

                if(!do_init && multi_atm_before != multi_atm_now) {
                  for(size_t lane = 0; lane < numLanes; ++lane) {
                    gpuLane_t* gpu = &(gpuLanes[lane]); // abbreviation
                    // Deallocation of GPU memory
                    free_memory_on_GPU((void **) &gpu->aero_G);
                    free_memory_on_GPU((void **) &gpu->obs_G);
                    free_memory_on_GPU((void **) &gpu->atm_G);
                    free_memory_on_GPU((void **) &gpu->tsurf_G);
                    free_memory_on_GPU((void **) &gpu->np_G);
                    free_memory_on_GPU((void **) &gpu->los_G);
                    free_memory_on_GPU((void **) &gpu->aero_beta_G);
                    if(multi_atm_before)
                      free_memory_on_GPU((void **) &gpu->atm_id_G);
                  }
                  free(gpuLanes);
                  debug_printf("[INFO] multi_atm_now = %d\n", (int) multi_atm_now);
                  debug_printf("[INFO] GPU lanes memory deallocation\n");
                }

                if(do_init) {
                  // Initialize ctl and tbl-struct (1 per GPU)
                  ctl_G = malloc_GPU(ctl_t, 1);
                  copy_data_to_GPU(ctl_G, ctl, sizeof(ctl_t), 0);
                  double tic = omp_get_wtime();
                  tbl_G = get_tbl_on_GPU(ctl);
                  double toc = omp_get_wtime();
                  printf("TIMER #%d jurassic-gpu reading table time: %lf\n",
                  ctl->MPIglobrank, toc - tic);
                }

                // Get number of possible lanes
                size_t gpuMemFree, gpuMemTotal;
                cuCheck(cudaMemGetInfo(&gpuMemFree, &gpuMemTotal));
                debug_printf("[INFO] memory GPU: free %.3f of total %.3f MByte = %.1f %%\n",
                              1e-6*gpuMemFree, 1e-6*gpuMemTotal, gpuMemFree/(.01*gpuMemTotal));

                numLanes = (size_t)((0.9*gpuMemFree) / (double)sizePerLane); // Only use 90% of free GPU memory ...
                printf("DEBUG #%d max possible number of Lanes: %d\n",
                ctl->MPIglobrank, numLanes);
                // ... other space is needed for alignment and profiling buffers
                size_t const maxNumLanes = 4; // Do not really need more than a handfull of lanes
                if (numLanes > maxNumLanes) numLanes = maxNumLanes;
                debug_printf("[INFO] GPU memory per lane: %.3f MByte, try to fit %i lanes\n", 1e-6*sizePerLane, numLanes);
                if (numLanes < 1) ERRMSG("Memory requirement per lane is too high, no lanes");

                gpuLanes = (gpuLane_t*) malloc(numLanes*sizeof(gpuLane_t)); // (this memory is never freed)
                for(size_t lane = 0; lane < numLanes; ++lane) {
                  gpuLane_t* gpu = &(gpuLanes[lane]); // abbreviation
                  // Allocation of GPU memory
                  gpu->aero_G		   = malloc_GPU(aero_t, 1);
                  gpu->obs_G		   = malloc_GPU(obs_t, 1);
                  gpu->tsurf_G	   = malloc_GPU(double, NR);
                  gpu->np_G		     = malloc_GPU(int, NR);
                  gpu->los_G		   = (pos_t (*)[NLOS])__allocate_on_GPU(NR*NLOS*sizeof(pos_t), __FILE__, __LINE__); 
                  gpu->aero_beta_G = (double(*)[ND])__allocate_on_GPU(NLMAX*ND*sizeof(double), __FILE__, __LINE__);
                  if(multi_atm_now) {
                    gpu->atm_id_G  = malloc_GPU(int, NR);
                    gpu->atm_G     = malloc_GPU(atm_t, NR);
                  }
                  else {
                    gpu->atm_id_G  = NULL;
                    gpu->atm_G	   = malloc_GPU(atm_t, 1);
                  }
                  cuCheck(cudaStreamCreate(&gpu->stream));
                  debug_printf("[INFO] cudaStreamCreate --> streamId %d\n", gpu->stream);
                } // lane
                debug_printf("[INFO] multi_atm_now = %d\n", (int) multi_atm_now);
                debug_printf("[INFO] GPU lanes memory allocation\n");
                cuCheck(cudaMemGetInfo(&gpuMemFree, &gpuMemTotal));
                debug_printf("[INFO] memory GPU: free %.3f of total %.3f MByte = %.1f %%\n",
                              1e-6*gpuMemFree, 1e-6*gpuMemTotal, gpuMemFree/(.01*gpuMemTotal));
              } // checkmode

				do_init = false;
        multi_atm_before = multi_atm_now;
        double toc = omp_get_wtime();
        printf("TIMER #%d jurassic-gpu gpu_lanes initialization time: %lf\n", ctl->MPIglobrank, toc - tic);
			} // do_init || (!do_init && multi_atm_before != multi_atm_now)
    } //<------- omp critical is here, maybe I should put it to the end of the function!

    // TODO: this may be a probelm in juwels-booster case (with more then one GPU  device per node)
    cudaSetDevice(0);

    if (ctl->checkmode) { printf("# %s: no operation in checkmode\n", __func__); return; }
    
    printf("numDevices: %d\n", numDevices);
    printf("DEBUG #%d numDevices: %d\n", ctl->MPIglobrank, numDevices);

    // I had add to this because of CPUs converting los to pos..    
    omp_set_nested(true);
#pragma omp parallel num_threads(numLanes)
#pragma omp for schedule(dynamic, 1) //work stealing
    for(int i = 0; i < n; i++) //loop over packages
    {
      int const myLane = omp_get_thread_num();
      assert(myLane < numLanes);
      char mask[NR][ND];
      save_mask(mask, &obs[i], ctl);
      copy_data_to_GPU(ctl_G, ctl, sizeof(ctl_t), gpuLanes[myLane].stream); // controls might change, update
      if(multi_atm_now)
        formod_one_package_GPU(ctl, ctl_G, tbl_G, divided_atms[i], &obs[i], aero, divided_atm_ids[i], &gpuLanes[myLane]);
      else
        formod_one_package_GPU(ctl, ctl_G, tbl_G, atm, &obs[i], aero, atm_id, &gpuLanes[myLane]);
      apply_mask(mask, &obs[i], ctl);
    }
    omp_set_nested(false);
    if(multi_atm_now) {
      for(int i = 0; i < n; i++) {
        free(divided_atms[i]);
        free(divided_atm_ids[i]);
      }
      free(divided_atms);
      free(divided_atm_ids);
    }
	} // jur_formod_multiple_packages_GPU
