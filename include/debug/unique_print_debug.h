//
// Created by Matteo Ranzi on 20/10/25.
//

#ifndef BM_P2P_COMM_PRINT_DEBUG_H
#define BM_P2P_COMM_PRINT_DEBUG_H

#include "print_debug.h"

#ifdef _OPENMP
#define PRAGMA_OMP_MASTER _Pragma("omp master")
#else
#define PRAGMA_OMP_MASTER
#endif

#define MASTER_PROCESS 0

#define UNIQUE_PRINT(rank, format, ...) do { if ((rank) == MASTER_PROCESS) PRAGMA_OMP_MASTER {printf("[UNIQUE] " format, ##__VA_ARGS__);} } while (0)

#define UNIQUE_PRINT_DEBUG_INFO(rank, format, ...) do { if ((rank) == MASTER_PROCESS) PRAGMA_OMP_MASTER {PRINT_DEBUG_INFO(rank, "[UNIQUE] " format, ##__VA_ARGS__);} } while (0)
#define UNIQUE_PRINT_DEBUG_WARN(rank, format, ...) do { if ((rank) == MASTER_PROCESS) PRAGMA_OMP_MASTER {PRINT_DEBUG_WARN(rank, "[UNIQUE] " format, ##__VA_ARGS__);} } while (0)
#define UNIQUE_PRINT_DEBUG_ERROR(rank, format, ...) do { if ((rank) == MASTER_PROCESS) PRAGMA_OMP_MASTER {PRINT_DEBUG_ERROR(rank, "[UNIQUE] " format, ##__VA_ARGS__);} } while (0)

#endif //BM_P2P_COMM_PRINT_DEBUG_H