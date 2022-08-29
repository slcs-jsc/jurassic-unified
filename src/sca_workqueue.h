#ifndef WORKQUEUE_H
#define WORKQUEUE_H

#include "jurassic.h"

/* Module containing all variables and functions related to a work-queue for scattering */

#define Queue_Inactive -1
#define Queue_Prepare   2
#define Queue_Execute   8
#define Queue_Collect  32

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

int jur_sca_init_queue(queue_t *q, int size);
int jur_sca_push_queue(queue_t *q, void* out, int ir);
int jur_sca_get_queue_item(queue_t *q, void **out, int *ir, int index);
int jur_sca_pop_queue(queue_t *q, void **out, int *ir);

#endif
