#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>

void dummyLinkSubproblem(int lId) // 修改 double dummyPenalty = 0.001;
{

  // *******1*****参数初始化
  CLink *pDummyLink;
  pDummyLink = &m_Link[lId]; // m_Link集合中包含 dummylink
  double temp_olf;
  double *update_OLPrimal_1st;
  double *update_OLPrimal_1st_cons;
  update_OLPrimal_1st = new double[m_Origin.size()];
  update_OLPrimal_1st_cons = new double[m_Origin.size()];

  double TempDummyLinkCost, TempDummyLinkCostDiff, snd_dummylink;
  double tempODdemand;
  double tempDummyLinkFlow;

  double o_Lp1std_linear;
  double tail_o_succecons = 0, head_o_succecons = 0;
  double tail_o_predecons = 0, head_o_predecons = 0;
  double tail_o_cons = 0, head_o_cons = 0, conscons = 0;

  // 2022.3.17 得到从基于destination的link flow 得到 OD demand，然后得到 dummy link flow
  // 基于origin的link子问题，ODdemand = sum( head(linkflow)=des)
  tempODdemand = getVariableODdemand(pDummyLink->pInNode->ID, pDummyLink->pOutNode->ID); // 大于等于0

  for (int oId = 0; oId < m_Origin.size(); oId++) // 对于每个起点
  {
    // *******2*****dummy link的一阶导，及目标函数的二阶导

    // 初始化origin based link flow，link flow

    tempDummyLinkFlow = fabs(m_ONUpperDem[oId][pDummyLink->pOutNode->ID]) - tempODdemand;

    if (tempDummyLinkFlow < 0) // 非负约束， 2022.03.27
    {
      tempDummyLinkFlow = 0;
    }

    // 更新dummy路段费用
    if (m_ONUpperDem[oId][pDummyLink->pOutNode->ID] == 0)
    {
      TempDummyLinkCost = InitDummyLinkFFT;
      TempDummyLinkCostDiff = 0;
    }
    else
    { // 更新dummy路段费用
      TempDummyLinkCost = InitDummyLinkFFT +
                          (1 / rstep) * (log(tempDummyLinkFlow / tempODdemand));

      TempDummyLinkCostDiff = (1 / rstep) * (tempODdemand / tempDummyLinkFlow); // error 分母为0
    }

    // 得到dummy路段的二阶导
    snd_dummylink = TempDummyLinkCostDiff + 2 * dummyPenalty;
    cout << "1. This link's 1st derivative (linkcost): " << TempDummyLinkCost << ". 2nd derivative: " << snd_dummylink << endl;

    // first oder derivative, link and origin based, second part, with oblf variable;
    o_Lp1std_linear = m_OLPrimal[oId][lId] * dummyPenalty * 2;

    // ******3.1******dummy link的pInNode节点的link残差
    tail_o_succecons = 0;
    for (int lId2 = 0; lId2 < pDummyLink->pInNode->OutgoingLink.size(); lId2++)
    {
      if (pDummyLink->pInNode->OutgoingLink[lId2] != lId) // 不包含 dummy link本身
      {
        tail_o_succecons += m_OLPrimal[oId][pDummyLink->pInNode->OutgoingLink[lId2]]; // 累加除待解决变量之外的origin based link flow（outgoing）
      }
    }
    tail_o_predecons = 0;
    for (int lId2 = 0; lId2 < pDummyLink->pInNode->IncomingLink.size(); lId2++)
      if (pDummyLink->pInNode->IncomingLink[lId2] != lId) // 不包含 dummy link本身
      {
        tail_o_predecons += m_OLPrimal[oId][pDummyLink->pInNode->IncomingLink[lId2]]; // 累加除待解决变量之外的origin based link flow（incoming）
      }
    tail_o_cons = tail_o_succecons - tail_o_predecons;

    tail_o_cons -= m_ONUpperDem[oId][pDummyLink->pOutNode->ID]; // pDummyLink->pInNode->ID表示destination，

    // cout<<" pLink_InNode_o_Outgoing - pLink_InNode_o_Incoming-ONDem: "<<tail_o_cons<<endl;
    // head node=out node

    // *******3.2*****dummy link的pOutNode节点的link残差
    head_o_succecons = 0;
    for (int lId2 = 0; lId2 < pDummyLink->pOutNode->OutgoingLink.size(); lId2++)
    {
      // cout<<" head_outgonging_links   "<<lId2<<", "<<pDummyLink->pOutNode->OutgoingLink[lId2]<<": ("<< m_Link[pDummyLink->pOutNode->OutgoingLink[lId2]].pInNode->Name <<", "<<m_Link[pLink->pOutNode->OutgoingLink[lId2]].pOutNode->Name<<")"<<endl;
      if (pDummyLink->pOutNode->OutgoingLink[lId2] != lId) // 不包含 dummy link本身
      {
        head_o_succecons += m_OLPrimal[oId][pDummyLink->pOutNode->OutgoingLink[lId2]];
        // cout<<" -> "<<m_OLPrimal[oId][pDummyLink->pOutNode->OutgoingLink[lId2]]<<endl;
      }
    }
    // cout<<" pLink_OutNode_o_Outgoing Link Flow Sum WO itself: "<<head_o_succecons<<endl;
    head_o_predecons = 0;
    for (int lId2 = 0; lId2 < pDummyLink->pOutNode->IncomingLink.size(); lId2++)
      if (pDummyLink->pOutNode->IncomingLink[lId2] != lId) // 不包含 dummy link本身
      {
        head_o_predecons += m_OLPrimal[oId][pDummyLink->pOutNode->IncomingLink[lId2]];
      }
    // cout<<" pLink_OutNode_o_Incoming Link Flow Sum WO itself: "<<head_o_predecons<<endl;
    head_o_cons = head_o_succecons - head_o_predecons;

    // 2022.3.11 将m_ONDem 改为 m_ONUpperDem
    head_o_cons -= m_ONUpperDem[oId][pDummyLink->pOutNode->ID]; // pDummyLink->pInNode->ID表示destination

    // cout<<" pLink_OutNode_o_Outgoing - pLink_OutNode_o_Incoming-ONDem: "<<head_o_cons<<endl;
    conscons = tail_o_cons - head_o_cons; // 对偶的常数部分

    // *******4***** 拟牛顿公式，更新origin-based link flow
    update_OLPrimal_1st_cons[oId] = dummyPenalty * conscons + m_ONDual[oId][pDummyLink->pInNode->ID] - m_ONDual[oId][pDummyLink->pOutNode->ID];
    // 目标函数的一阶导 d = link_BRP + 2*penalty*(v_o_a)+(r_tail-r_head)+penalty*(e_tail-e_head)

    double diffDual = m_ONDual[oId][pDummyLink->pInNode->ID] - m_ONDual[oId][pDummyLink->pOutNode->ID];

    // 目标函数的一阶导
    update_OLPrimal_1st[oId] = TempDummyLinkCost + o_Lp1std_linear + update_OLPrimal_1st_cons[oId];

    // 更新第j次基于起点的路段流量

    m_OLPrimal[oId][lId] = m_OLPrimal[oId][lId] - (dummy_step * update_OLPrimal_1st[oId]) / snd_dummylink;
    // m_OLPrimal一直等于0， d/s一直大于0，由于非负约束，update m_OLPrimal等于0

    if (m_OLPrimal[oId][lId] < 0)
      m_OLPrimal[oId][lId] = 0;

    printf("2  m_OLPrimal[%d][%d]:%f ,update_OLPrimal_1st[oId]:%f,snd_dummylink:%f\n",
           oId, lId, m_OLPrimal[oId][lId], update_OLPrimal_1st[oId], snd_dummylink);
  }

  LinkFlow[lId] = tempDummyLinkFlow; // 2022.03.23

  double terminate, snd_15, d_15, terminategap = 10, terminate_0 = 10; // termination criteria
  int counter = 1;
  terminate = 0;
  for (int oId = 0; oId < m_Origin.size(); oId++)
    terminate += fabs(m_OLPrimal[oId][lId] * update_OLPrimal_1st[oId]); // flow(k+1)*gradient(k)
                                                                        // terminate+=fabs(update_OLPrimal[oId]*update_OLPrimal_1st[oId]);

  // dummylink的term条件需要单独设置
  int flag = 0;
  while (flag == 1 && terminate > 1e-5 && terminategap > 1e-3) // terminate > 1e-9   //counter<1  //
  // 这里用了两个指标
  {
    counter = counter + 1;
    for (int oId = 0; oId < m_Origin.size(); oId++)
    {
      // *******5*****再次更新 dummy link的一阶导，及目标函数的二阶导
      // 更新DummyLinkFlow
      tempODdemand = getVariableODdemand(pDummyLink->pInNode->ID, pDummyLink->pOutNode->ID); // 大于等于0
      // printf("*****2*** dummylink:%d, tempODdemand:%f \n", lId, tempODdemand);
      tempDummyLinkFlow = fabs(m_ONUpperDem[oId][pDummyLink->pOutNode->ID]) - tempODdemand;
      if (tempDummyLinkFlow < 0) // 非负约束， 2022.03.27
      {
        tempDummyLinkFlow = 0;
      }

      if (m_ONUpperDem[oId][pDummyLink->pOutNode->ID] == 0)
      {
        TempDummyLinkCost = InitDummyLinkFFT;
        TempDummyLinkCostDiff = 0;
      }
      else
      {
        // 更新dummy路段费用
        TempDummyLinkCost = InitDummyLinkFFT +
                            (1 / rstep) * (log(tempDummyLinkFlow / tempODdemand));

        TempDummyLinkCostDiff = (1 / rstep) * (tempODdemand / tempDummyLinkFlow);
      }

      snd_dummylink = TempDummyLinkCostDiff + 2 * dummyPenalty; // 由于link flow是变化的，snd也会变化，

      update_OLPrimal_1st[oId] = TempDummyLinkCost + m_OLPrimal[oId][lId] * dummyPenalty * 2 + update_OLPrimal_1st_cons[oId]; // origin-based link flow会变化，常数项不变了；
                                                                                                                              // 更新第j+1次基于起点的路段流量
      // *******6***** 再次拟牛顿公式，更新origin-based link flow
      temp_olf = m_OLPrimal[oId][lId] - (dummy_step * update_OLPrimal_1st[oId]) / snd_dummylink; /// 马上update，接着做投影；
      if (temp_olf < 0)
        temp_olf = 0;

      tempDummyLinkFlow -= m_OLPrimal[oId][lId]; // 减去第j次基于起点的路段流量
      tempDummyLinkFlow += temp_olf;             // 这里只更新一个，很聪明的做法；//加上第j+1次 更新 基于起点的路段流量

      m_OLPrimal[oId][lId] = temp_olf;

      printf(" innerloop:%d,  DummyLinkFlow:%f, m_OLPrimal[1][6]:%f, terminate:%f, terminategap:%f \n",
             counter, tempDummyLinkFlow, m_OLPrimal[oId][lId], terminate, terminategap);
    }

    terminate = 0;
    for (int oId = 0; oId < m_Origin.size(); oId++)
      terminate += fabs(update_OLPrimal_1st[oId] * m_OLPrimal[oId][lId]);
    // cout<<"  link based terminate: "<<terminate<<endl;
    terminategap = fabs(terminate - terminate_0);
    terminate_0 = terminate;
  }
  LinkFlow[lId] = tempDummyLinkFlow;
  delete[] update_OLPrimal_1st;
  delete[] update_OLPrimal_1st_cons;
  cout << "dummy AG:" << setprecision(20) << terminate << endl;

  cout << "      dummyLink in Block: (" << m_Link[lId].pInNode->Name << "->" << m_Link[lId].pOutNode->Name
       << "), dummyflow=" << tempDummyLinkFlow << ",m_OLPrimal[1][6]:" << m_OLPrimal[0][lId] << endl;
  cout << endl;
}
