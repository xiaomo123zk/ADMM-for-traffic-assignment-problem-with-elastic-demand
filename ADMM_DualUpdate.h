void ADMM_DualUpdate()
{
  printf("\n\n***********ADMM_DualUpdate\n");
  double oId_nId_outgoingLfCons = 0;
  double oId_nId_IncomingLfCons = 0;
  double on_conscons = 0;
  int link_tail;
  CNode *pNode;
  CLink *pLink;

  for (int oId = 0; oId < m_Origin.size(); oId++) // 对于每个origin
  {
    // m_Node.size() 替换为2022.3.13 m_Node_withDummyLink.size()
    for (int nId = 0; nId < m_Node_withDummyLink.size(); nId++) // 遍历每个节点的 所有OutgoingLink和IncomingLink
    {
      // pNode = &m_Node.at(nId);
      pNode = &m_Node_withDummyLink.at(nId);

      // 关于节点流量守恒时，以Small net为例，应该考虑dummy link，2022-03，怎么添加
      // 构建node与dummy link的关系，已知dummy link只与origin node和destination node有关
      oId_nId_outgoingLfCons = 0;
      if (pNode->OutgoingLink.size() > 0)
        for (int lId = 0; lId < pNode->OutgoingLink.size(); lId++) // 遍历该节点的所有OutgoingLinks
        {
          oId_nId_outgoingLfCons += m_OLPrimal[oId][pNode->OutgoingLink[lId]];
          // printf("Outlink   [oId:%d][nId:%d], m_OLPrimal[%d][Outlink:%d]:%f \n",
          //        oId + 1, nId + 1, oId + 1, pNode->OutgoingLink[lId] + 1, m_OLPrimal[oId][pNode->OutgoingLink[lId]]);
        }

      oId_nId_IncomingLfCons = 0;
      if (pNode->IncomingLink.size() > 0)
        for (int lId = 0; lId < pNode->IncomingLink.size(); lId++) // 遍历该节点的所有IncomingLinks
        {
          oId_nId_IncomingLfCons += m_OLPrimal[oId][pNode->IncomingLink[lId]];
          // printf("Inlink   [oId:%d][nId:%d], m_OLPrimal[%d][Inlink:%d]:%f \n",
          //        oId + 1, nId + 1, oId + 1, pNode->IncomingLink[lId] + 1, m_OLPrimal[oId][pNode->IncomingLink[lId]]);
        }

      // on_conscons = oId_nId_outgoingLfCons - oId_nId_IncomingLfCons - m_ONDem[oId][nId];
      // 残差 = sum(节点的outgoingLink) - sum(节点的IncomingLink)- g_o_d
      on_conscons = oId_nId_outgoingLfCons - oId_nId_IncomingLfCons - m_ONUpperDem[oId][nId];

      m_ONDual[oId][nId] += penalty * on_conscons; // 不断累加的过程
    }
  }
}
