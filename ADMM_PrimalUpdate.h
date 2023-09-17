#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>

void ADMM_PrimalUpdate()
{

  // s_primal = omp_get_wtime();
  CLink *pLink;
  for (int bId = 0; bId <= max_BlockId; bId++)
  {
    // cout<<"BlockId:  "<<bId+1<<endl;
    max_period = 0;
    double tempODdemand = 0;

    int lId;
#pragma omp parallel private(lId) num_threads(num_threads2)
    {
#pragma omp for
      for (lId = 0; lId < BlockSet[bId].size(); lId++)
      {
        // s_primal_eff = omp_get_wtime();
        s_primal_eff = clock();

        int tempLink = BlockSet[bId][lId];

        pLink = &m_Link[tempLink];

        if (pLink->DummyLinkFlag == 1) // 是dummy link
        {
          printf("\n\n*************** dummylink start\n");
          dummyLinkSubproblem(BlockSet[bId][lId]); // 10
          printf("*************** dummylink end\n");
        }
        else // 不是 dummy link
        {
          printf("\n\n*****auto*start*********  block:%d's link:%d \n", bId + 1, BlockSet[bId][lId] + 1);
          linkSubproblem_v3(BlockSet[bId][lId]); // 10
          printf("*****auto*end*********  block:%d's link:%d \n", bId + 1, BlockSet[bId][lId] + 1);
        }
        e_primal_eff = omp_get_wtime();

        if (max_period < (e_primal_eff - s_primal_eff))
          max_period = (e_primal_eff - s_primal_eff);
      }
    }
    p_primal_eff += max_period;
  }
}