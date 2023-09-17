#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
#define INF 1.0e13
#define ZERO 1.0e-13
#define INVALID -1      /* Represents an invalid value. */
#define WAS_IN_QUEUE -7 /* Shows that the node was in the queue before. (7 is for luck.) */

using namespace std;
class CNode // 节点
{
public:
  int ID;
  string Name;
  int Origin_ID = -1;
  vector<int> IncomingLink;
  vector<int> OutgoingLink;
};
class CLink // LinkID --> Fft，Capacity, BlockID, Alpha, Power
{
public:
  int ID;
  CNode *pInNode;
  CNode *pOutNode;
  double FreeFlowTravelTime;
  double Capacity;
  double Alpha = 0.15;
  double Power = 4.0;
  int BlockID;
  int DummyLinkFlag = 0;
};
class CODPair // ODPairID-->
{
public:
  int ori, des, desID, oriID;
  double ODdemand;
};
class COrigin
{
public:
  int ID;
  CNode *pOriginNode;
  vector<int> DestinationNode;
  vector<double> ODDemand;
  vector<int> ODPairIndex;
};
