double GetlinkUEGapED(int counter)
{
  CLink *pLink;
  double link_systemcost = 0.0, od_systemcost = 0.0;
  int *temp = new int[1]{6};
  int id, ori = 0;

  double sum_dummlylink = 0;
  double tempODdemand, tempDummyLinkFlow, TempDummyLinkCost;
  double diff_link_systemcost = 0.0, diff_linkcost = 0.0, diff_linkflow = 0.0;
  double dummylink_systemcost = 0.0, diff_dummylink_systemcost = 0.0, diff_dummylinkcost = 0.0, diff_dummylinkflow = 0.0;

  // ***********fenmu**************
  for (int link = 0; link < m_Link.size(); link++)
  {
    pLink = &m_Link[link];
    if (pLink->DummyLinkFlag != 1)
      link_systemcost += LinkFlow[link] * LinkCost[link];
  }

  for (int dummy = 0; dummy < 1; dummy++)
  {
    id = temp[dummy];
    pLink = &m_Link[id - 1];

    if (pLink->DummyLinkFlag == 1)
    {
      tempODdemand = getVariableODdemand(pLink->pInNode->ID, pLink->pOutNode->ID);
      tempDummyLinkFlow = fabs(m_ONUpperDem[ori][pLink->pOutNode->ID]) - tempODdemand;
      if (tempDummyLinkFlow < 0)
        tempDummyLinkFlow = 0;

      if (m_ONUpperDem[ori][pLink->pOutNode->ID] == 0)
      {
        TempDummyLinkCost = InitDummyLinkFFT;
      }
      else
      { // 更新dummy路段费用

        TempDummyLinkCost = InitDummyLinkFFT +
                            (1 / rstep) * (log(tempDummyLinkFlow / tempODdemand));
      }

      sum_dummlylink += TempDummyLinkCost * tempDummyLinkFlow;

      if (counter % 2 == 1) // 前后两次 当nOutloop是奇数时，就进入
      {

        DummyOddLinkCost_1_4 = TempDummyLinkCost;
        DummyOddLinkFlow_1_4 = tempDummyLinkFlow;
      }
      else if (counter % 2 == 0 || counter == 1)
      {

        DummyEvenLinkCost_1_4 = TempDummyLinkCost;
        DummyEvenLinkFlow_1_4 = tempDummyLinkFlow;
      }

      // diff_dummylink_systemcost = DummyEvenLinkFlow_1_4 * fabs(DummyEvenLinkCost_1_4 - DummyOddLinkCost_1_184) ;
      diff_dummylink_systemcost = DummyEvenLinkCost_1_4 * fabs(DummyEvenLinkFlow_1_4 - DummyOddLinkFlow_1_4);
    }
  }

  // ***********fenzi**************
  for (int link = 0; link < m_Link.size(); link++)
  {
    pLink = &m_Link[link];
    if (pLink->DummyLinkFlag != 1)
    {
      diff_linkflow = fabs(previousEvenLinkFlow[link] - previousOddLinkFlow[link]);
      diff_link_systemcost += LinkCost[link] * diff_linkflow;
    }
  }

  double fenzi = diff_link_systemcost + diff_dummylink_systemcost;
  double fenmu = link_systemcost + sum_dummlylink;
  delete[] temp;
  return fenzi / fenmu;
}
