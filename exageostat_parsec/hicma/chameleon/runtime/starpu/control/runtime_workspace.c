/**
 *
 * @file runtime_workspace.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU workspaces routines
 *
 * @version 1.0.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @date 2011-06-01
 *
 */
#include "chameleon_starpu.h"

static void RUNTIME_allocate_workspace_on_workers(void *arg)
{
    struct morse_starpu_ws_s *workspace = arg;
    enum starpu_worker_archtype type = 0;
    (void)type;

    int id = starpu_worker_get_id();

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
    type = starpu_worker_get_type(id);
    if (type == STARPU_CUDA_WORKER)
    {
        int memory_location = workspace->memory_location;

        if (memory_location == MORSE_HOST_MEM)
        {
            /* Use pinned memory because the kernel is very likely
             * to transfer these data between the CPU and the GPU.
             * */
            cudaMallocHost(&workspace->workspaces[id], workspace->size);
        }
        else {
            /* Allocate on the device */
            cudaMalloc(&workspace->workspaces[id], workspace->size);
        }
    }
    else
#endif
    {
        /* This buffer should only be used within the CPU kernel, so
         * there is no point in using pinned memory here. */
        workspace->workspaces[id] = malloc(workspace->size);
    }

    assert(workspace->workspaces[id]);
}


static void RUNTIME_free_workspace_on_workers(void *arg)
{
    struct morse_starpu_ws_s *workspace = arg;
    enum starpu_worker_archtype type = 0;
    (void)type;
    int id = starpu_worker_get_id();

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
    type = starpu_worker_get_type(id);
    if (type == STARPU_CUDA_WORKER)
    {
        int memory_location = workspace->memory_location;

        if (memory_location == MORSE_HOST_MEM)
        {
            cudaFreeHost(workspace->workspaces[id]);
        }
        else {
            cudaFree(workspace->workspaces[id]);
        }
    }
    else
#endif
    {
        free(workspace->workspaces[id]);
    }

    workspace->workspaces[id] = NULL;
}

/*
 * This function creates a workspace on each type of worker in "which_workers"
 * (eg. MORSE_CUDA|MORSE_CPU for all CPU and GPU workers).  The
 * memory_location argument indicates whether this should be a buffer in host
 * memory or in GPU memory (MORSE_HOST_MEM or MORSE_GPU_MEM). This function
 * returns 0 upon successful completion.:
 */
int RUNTIME_starpu_ws_alloc(MORSE_starpu_ws_t **workspace,
                            size_t size, int which_workers, int memory_location)
{
    if (!workspace)
        return -EINVAL;

    struct morse_starpu_ws_s *descr = calloc(1, sizeof(struct morse_starpu_ws_s));

    *workspace = descr;

    if (!descr)
        return -ENOMEM;

    descr->size = size;
    descr->memory_location = memory_location;
    descr->which_workers = which_workers;

    starpu_execute_on_each_worker(RUNTIME_allocate_workspace_on_workers, descr, which_workers);

    return 0;
}

int RUNTIME_starpu_ws_free(MORSE_starpu_ws_t *workspace)
{
    if (!workspace)
        return -EINVAL;

    starpu_execute_on_each_worker(RUNTIME_free_workspace_on_workers, workspace, workspace->which_workers);

    free(workspace);

    return 0;
}

void *RUNTIME_starpu_ws_getlocal(MORSE_starpu_ws_t *workspace)
{
    struct morse_starpu_ws_s *descr = workspace;
    int id = starpu_worker_get_id();
    return descr->workspaces[id];
}
