#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<sys/time.h>
#include<omp.h>
#include"utils.h"
#include"input.h"
#include"mesh.h"
#include"solver.h"
#include"flux.h"
#include"boundary.h"
#include"sa.h"

int main(int argc, char **argv)
{

    struct timeval start, stop;             // 用于计算运行时间
    
    gettimeofday(&start,NULL);              // 得到当前时间，秒和微秒     
   
    SOLVER* solver = solverInit(argv[1]);   // 初始化Solver
    
    solverSolve(solver);                    // 
    
    solverWriteSolution(solver);

    solverWriteRestart(solver);
      
    solverFree(solver); 

    gettimeofday(&stop,NULL);
    
    printf("\nDuration %f s\n", duration(start, stop));

    return 0;

}
