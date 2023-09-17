#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>

void linkSubproblem(int lId)
{

  CLink *pLink;

  pLink = &m_Link[lId]; //???
  double temp_olf;
  double *update_OLPrimal_1st;
  double *update_OLPrimal_1st_cons;
  update_OLPrimal_1st = new double[m_Origin.size()];
  update_OLPrimal_1st_cons = new double[m_Origin.size()];
  // 初始化origin based link flow，link flow
  double templinkFlow = LinkFlow[lId];

  cout << "0. Link in Block: (" << m_Link[lId].pInNode->Name << "->" << m_Link[lId].pOutNode->Name << "), flow=" << templinkFlow << endl;
  // 计算2nd，1st 导数
  double temppow = templinkFlow * templinkFlow * templinkFlow;
  double snd = LinkCostDiffCoef[lId] * temppow + 2 * penalty;
  if (lId == 3)
    cout << ",link flow:" << templinkFlow << "1. This link's 2st derivative (linkcost): " << LinkCostDiffCoef[lId] * temppow << endl;

  LinkCost[lId] = LinkFreeTravelTime[lId] + LinkCostCoef[lId] * temppow * templinkFlow;
  if (lId == 3)
    cout << "1. This link's 1st derivative (linkcost): " << LinkCost[lId] << ". 2nd derivative: " << snd << ", linkcostdiff: " << LinkCostDiffCoef[lId] << endl;
  double o_Lp1std_linear;
  double tail_o_succecons = 0, head_o_succecons = 0;
  double tail_o_predecons = 0, head_o_predecons = 0;
  double tail_o_cons = 0, head_o_cons = 0, conscons = 0;

  for (int oId = 0; oId < m_Origin.size(); oId++)
  {
    // first oder derivative, link and origin based, second part, with oblf variable;
    o_Lp1std_linear = m_OLPrimal[oId][lId] * penalty * 2;
    tail_o_succecons = 0;
    for (int lId2 = 0; lId2 < pLink->pInNode->OutgoingLink.size(); lId2++)
      if (pLink->pInNode->OutgoingLink[lId2] != lId)                             //
        tail_o_succecons += m_OLPrimal[oId][pLink->pInNode->OutgoingLink[lId2]]; // 包含dummy link，累加除待解决变量之外的origin based link flow（outgoing）
    tail_o_predecons = 0;
    for (int lId2 = 0; lId2 < pLink->pInNode->IncomingLink.size(); lId2++)
      if (pLink->pInNode->IncomingLink[lId2] != lId)
        tail_o_predecons += m_OLPrimal[oId][pLink->pInNode->IncomingLink[lId2]]; // 包含dummy link，累加除待解决变量之外的origin based link flow（incoming）
    tail_o_cons = tail_o_succecons - tail_o_predecons;

    tail_o_cons -= m_ONUpperDem[oId][pLink->pOutNode->ID]; // ?? 2022.03 m_ONUpperDem or m_ONDem

    // head node=out node
    head_o_succecons = 0;
    for (int lId2 = 0; lId2 < pLink->pOutNode->OutgoingLink.size(); lId2++)
    {
      if (pLink->pOutNode->OutgoingLink[lId2] != lId)
      {
        head_o_succecons += m_OLPrimal[oId][pLink->pOutNode->OutgoingLink[lId2]]; // 包含dummy link
      }
    }
    head_o_predecons = 0;
    for (int lId2 = 0; lId2 < pLink->pOutNode->IncomingLink.size(); lId2++)
      if (pLink->pOutNode->IncomingLink[lId2] != lId)
        head_o_predecons += m_OLPrimal[oId][pLink->pOutNode->IncomingLink[lId2]]; // 包含dummy link，
    head_o_cons = head_o_succecons - head_o_predecons;

    head_o_cons -= m_ONUpperDem[oId][pLink->pOutNode->ID]; // ?? 2022.03 m_ONUpperDem m_ONDem
    // cout << " pLink_OutNode_o_Outgoing - pLink_OutNode_o_Incoming-ONDem: " << head_o_cons << endl;
    conscons = tail_o_cons - head_o_cons; // 常数部分

    update_OLPrimal_1st_cons[oId] = penalty * conscons + m_ONDual[oId][pLink->pInNode->ID] - m_ONDual[oId][pLink->pOutNode->ID];

    double diffDual = m_ONDual[oId][pLink->pInNode->ID] - m_ONDual[oId][pLink->pOutNode->ID];
    update_OLPrimal_1st[oId] = LinkCost[lId] + o_Lp1std_linear + update_OLPrimal_1st_cons[oId];

    m_OLPrimal[oId][lId] = m_OLPrimal[oId][lId] - (update_OLPrimal_1st[oId]) / snd;

    if (m_OLPrimal[oId][lId] < 0)
      m_OLPrimal[oId][lId] = 0;
  }
  // printf("-**-11- ,LinkFlow[%d]:%f, m_OLPrimal[0][%d]:%f \n", lId, LinkFlow[lId], lId, m_OLPrimal[0][lId]);

  double terminate, snd_15, d_15, terminategap = 10, terminate_0 = 10; // termination criteria
  int counter = 1;
  terminate = 0;
  for (int oId = 0; oId < m_Origin.size(); oId++)
    terminate += fabs(m_OLPrimal[oId][lId] * update_OLPrimal_1st[oId]); // flow(k+1)*gradient(k)
                                                                        // terminate+=fabs(update_OLPrimal[oId]*update_OLPrimal_1st[oId]);
  // cout <<"3. Loop: "<<endl;
  templinkFlow = 0;
  for (int oId = 0; oId < m_Origin.size(); oId++)
  {
    templinkFlow += m_OLPrimal[oId][lId];
  }

  // 每个origin都update了，所以所有的origin都需要重新相加；
  while (terminate > 1e-5 && terminategap > 1e-3) // terminate > 1e-9   //counter<1  //
  // 这里用了两个指标
  {
    // printf("innerloop:%d, terminate:%f, terminategap:%f \n", counter, terminate, terminategap);
    counter = counter + 1;
    for (int oId = 0; oId < m_Origin.size(); oId++)
    {
      temppow = templinkFlow * templinkFlow * templinkFlow;                                                          // 这里可以改成link flow[lid]
      snd = LinkCostDiffCoef[lId] * temppow + 2 * penalty;                                                           // 由于link flow是变化的，snd也会变化，
      LinkCost[lId] = LinkFreeTravelTime[lId] + LinkCostCoef[lId] * temppow * templinkFlow;                          // 1st中的BPR time也会随着link flow的变化而变化。
      update_OLPrimal_1st[oId] = LinkCost[lId] + m_OLPrimal[oId][lId] * penalty * 2 + update_OLPrimal_1st_cons[oId]; // origin-based link flow会变化，常数项不变了；

      // temp_olf=update_OLPrimal[oId]-update_OLPrimal_1st[oId]/snd;///
      temp_olf = m_OLPrimal[oId][lId] - update_OLPrimal_1st[oId] / snd; /// 马上update，接着做投影；
      if (temp_olf < 0)
        temp_olf = 0;

      templinkFlow -= m_OLPrimal[oId][lId];
      templinkFlow += temp_olf; // 这里只更新一个，很聪明的做法；
      // update_OLPrimal[oId]=temp_olf;
      m_OLPrimal[oId][lId] = temp_olf;
      // printf("-02-  after0  templinkFlow:%f, m_OLPrimal:%f\n",
      // templinkFlow, m_OLPrimal[oId][lId]);
    }
    terminate = 0;
    for (int oId = 0; oId < m_Origin.size(); oId++)
      terminate += fabs(update_OLPrimal_1st[oId] * m_OLPrimal[oId][lId]);
    // cout<<"  link based terminate: "<<terminate<<endl;
    terminategap = fabs(terminate - terminate_0);
    terminate_0 = terminate;
  }
  LinkFlow[lId] = templinkFlow;
  delete[] update_OLPrimal_1st;
  delete[] update_OLPrimal_1st_cons;

  cout << "      Link in Block: (" << m_Link[lId].pInNode->Name << "->" << m_Link[lId].pOutNode->Name << "), flow=" << templinkFlow << ",m_OLPrimal:" << m_OLPrimal[0][lId] << endl;
}
