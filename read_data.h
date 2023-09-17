#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
void ReadData(string fileName)
{
  int count;
  string lineStr;   // 读取每行的字符串
  string inlineStr; // 每行断开
  stringstream stream;
  int intData;
  double doubleData;
  // read node_set
  ifstream nodeFile(fileName + "_node.csv");
  if (nodeFile.is_open())
  {
    count = 0;
    while (getline(nodeFile, lineStr)) // 当没有读取到文件末尾时循环继续
    {
      CNode Node;
      stringstream ss(lineStr);
      getline(ss, inlineStr, '\t');
      Node.Name = inlineStr; // Node的name是连续的，从1开始，ID从0开始，间隔-1；
      Node.ID = count;
      count++;
      m_Node.push_back(Node);
      m_Node_withDummyLink.push_back(Node);
    }
    nodeFile.close();
    printf("Read node data finish!\n");
  }
  else
  {
    printf("Read node data wrong!\n");
  }

  // read link_set
  ifstream linkFile(fileName + "_link.csv"); //
  count = 0;
  if (linkFile.is_open())
  {
    while (getline(linkFile, lineStr))
    {
      CLink Link;
      stringstream ss(lineStr);
      Link.ID = count;
      count++;

      getline(ss, inlineStr, ',');
      stream << inlineStr;
      stream >> intData;
      m_Node[intData - 1].OutgoingLink.push_back(Link.ID);
      m_Node_withDummyLink[intData - 1].OutgoingLink.push_back(Link.ID); // 2022.03
      Link.pInNode = &m_Node[intData - 1];
      stream.clear();

      getline(ss, inlineStr, ',');
      stream << inlineStr;
      stream >> intData;
      m_Node[intData - 1].IncomingLink.push_back(Link.ID);
      m_Node_withDummyLink[intData - 1].IncomingLink.push_back(Link.ID); // 2022.03
      Link.pOutNode = &m_Node[intData - 1];
      stream.clear();

      getline(ss, inlineStr, ',');
      stream << inlineStr;
      stream >> doubleData;
      Link.FreeFlowTravelTime = doubleData;
      stream.clear();

      getline(ss, inlineStr, ',');
      stream << inlineStr;
      stream >> doubleData;
      Link.Capacity = doubleData;
      stream.clear();

      getline(ss, inlineStr, ',');
      stream << inlineStr;
      stream >> intData;
      Link.BlockID = intData - 1; // 要保证Block Name从1开始依次累加，与Block ID有固定的-1关系
      stream.clear();

      m_Link.push_back(Link);
    }
    linkFile.close();
    printf("Read link data finish!\n");
  }
  else
  {
    printf("Read link data wrong!\n");
  }

  // establish Block_set
  count = 0;
  max_BlockId = 0;
  for (int i = 0; i < m_Link.size(); i++)
    if (m_Link[i].BlockID > max_BlockId)
      max_BlockId = m_Link[i].BlockID;
  // cout<<max_BlockId<<endl;
  BlockSet = new vector<int>[max_BlockId + 1]; // 要保证Block Name与Block ID有固定的-1关系
  for (int i = 0; i < max_BlockId + 1; i++)
  {
    vector<int> LinkIdSet;
    for (int j = 0; j < m_Link.size(); j++)
      if (i == m_Link[j].BlockID)
        LinkIdSet.push_back(j);
    BlockSet[i] = LinkIdSet;
  }

  // read OD demand and origin set//////////////////////////////////
  // read fileName_od.csv
  ifstream odFile(fileName + "_od.csv");
  CNode *pNode;
  COrigin *pOrigin;
  int line = -1;
  if (odFile.is_open())
  {
    while (getline(odFile, lineStr)) // 当没有读取到文件末尾时循环继续
    {
      line++;
      stringstream ss(lineStr);
      getline(ss, inlineStr, ',');
      stream << inlineStr;
      stream >> intData;
      pNode = &m_Node[intData - 1]; // 要保证nodeName与nodeID有固定的-1关系
      stream.clear();

      if (pNode->Origin_ID == -1)
      {
        // pOrigin = new COrigin();
        COrigin Origin;
        Origin.ID = m_Origin.size(); // 也是要保持originID，这个很有问题。
        Origin.pOriginNode = pNode;
        pNode->Origin_ID = Origin.ID; // Origin.ID表示m_Origin中的Index，这个Index记录在了pNode中。
        m_Origin.push_back(Origin);   // node->m_ori  m_ori->node都有可以获取，双向键值对。
        pOrigin = &m_Origin[pNode->Origin_ID];
      }
      else
      {
        pOrigin = &m_Origin[pNode->Origin_ID];
      }
      getline(ss, inlineStr, ',');
      stream << inlineStr;
      stream >> intData;
      pNode = &m_Node[intData - 1]; // 这里获取的是des的信息，也是Node的ID与Name有固定关系。
      stream.clear();

      getline(ss, inlineStr, ',');
      stream << inlineStr;
      stream >> doubleData;

      if (doubleData > 0)
      {
        pOrigin->DestinationNode.push_back(pNode->ID); // 这里只是获得了Node的ID信息。
        pOrigin->ODDemand.push_back(doubleData);
        pOrigin->ODPairIndex.push_back(line); // 这里也是为了做双向索引。 ???
      }
      stream.clear();
    }
    odFile.close();
    printf("Read origin data finish!\n");
  }
  else
  {
    printf("Read origin data wrong!\n");
  }

  // Read O_N based Demand table
  m_ONUpperDem = new double *[m_Origin.size()]; // 2022.3.11 OD upperbound demand
  m_ONDem = new double *[m_Origin.size()];
  for (int oId = 0; oId < m_Origin.size(); oId++)
  {
    m_ONUpperDem[oId] = new double[m_Node.size()];
    m_ONDem[oId] = new double[m_Node.size()];
  }
  // ifstream onFile(fileName + "_allON_with_dummmy1-24.csv");
  ifstream onFile(fileName + "_on.csv");
  int oriId, desId;
  if (onFile.is_open())
  {
    while (getline(onFile, lineStr))
    {
      stringstream ss(lineStr);

      getline(ss, inlineStr, ',');
      stream << inlineStr;
      stream >> intData;
      oriId = intData - 1; // 获取的也是NodeID，前提是默认Node ID与Name相差-1
      stream.clear();

      getline(ss, inlineStr, ',');
      stream << inlineStr;
      stream >> intData;
      desId = intData - 1;
      stream.clear();

      getline(ss, inlineStr, ',');
      stream << inlineStr;
      stream >> doubleData;
      m_ONDem[oriId][desId] = doubleData;
      stream.clear();

      getline(ss, inlineStr, ',');
      stream << inlineStr;
      stream >> doubleData;
      m_ONUpperDem[oriId][desId] = doubleData;

      stream.clear();
      // cout << oriId << "," << desId << "," << m_ONDem[oriId][desId] << endl;
    }
    onFile.close();
    printf("Read on data finish!\n");
    // for(int oi=0;oi<m_Origin.size();oi++)
    //     for (int ni=0;ni<m_Node.size();ni++)
    //         cout<<oi+1<<","<<ni+1<<", "<<m_ONDem[oi][ni]<<endl;
  }
  else
  {
    printf("Read on data wrong!\n");
  }
}
