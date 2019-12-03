#include "hw5.h"
#include <map>
#include<iostream>

using namespace std;


int main()
{
  // initialize problem
  map<int, region> regions ;
  regions[1].setRegionID(1);
  regions[1].setScatter(45);
  regions[1].setAbs(4);
  regions[1].setLims({0,5});
  regions[2].setRegionID(2);
  regions[2].setScatter(8);
  regions[2].setAbs(50);
  regions[2].setLims({5,10});
  randNum r1;

  map<int,int> numAbs = {{1,0},{2,0}};
  map<int,int> numFScats ={{1,0},{2,0}};

  for (int iEv =0; iEv<100000; iEv++)
  {
    if (iEv%10000==0) cout << "current event "<<iEv<<"\n";
    /// initialize event   
    double init_loc = 5.*r1.next();
    while (init_loc==5) init_loc = 5.*r1.next();
    double init_dirn = 2.*r1.next()-1;

    particle part;
    part.setLoc(init_loc);
    part.setRegion(1);
    part.setDir(init_dirn);

    while (part.getStatus())
    {

      double next_loc = regions[part.getRegion()].getNextLoc(r1.next(), (&part));
      part.setLoc(next_loc);
      int currentRegion = regions[part.getRegion()].checkRegion((&part));
      if (currentRegion==part.getRegion())
      {
        int coll = regions[part.getRegion()].getColl(r1.next());
        if (coll==0) 
        {
          part.setStatus(false);
          numAbs[currentRegion]++;
        }
        else
        {
          double nextDir = 2.*r1.next()-1;
          while (nextDir==0) nextDir = 2.*r1.next()-1;
          part.setDir(nextDir);
          if (nextDir>0) numFScats[currentRegion]++;
        }
      }
      else 
      {
        std::map<int,region>::iterator it = regions.find(currentRegion);
        if (it==regions.end())
        {
          part.setStatus(false);
        }
        else
        {
          if (currentRegion<part.getRegion())
          {
            part.setLoc(regions[currentRegion].getLims()[1]);
          }
          else part.setLoc(regions[currentRegion].getLims()[0]);
          part.setRegion(currentRegion);
        }
      }
    }
  }

  cout << "Num Absorptions in region1 = " << numAbs[1] << "\n";
  cout << "Num Absorptions in region2 = " << numAbs[2] << "\n";
  cout << "Num forward scatters in region1 = " << numFScats[1] <<"\n";
  cout << "Num forward scatters in region2 = " << numFScats[2] <<"\n";

  return 0;
}