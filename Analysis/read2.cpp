// g++ -std=c++11 ./read.cpp -o ./read;
#include <iostream>
#include <fstream>
//#include<string.h>
using namespace std;

int main (int argc, char* argv[]) {
  string line, pre="./Au_p0.1_i17500_", end=".dump", str;
  ofstream output, fliquid;
  ifstream input, fentro;
  int multi = 1, start = 41;
  int id, nodesid, stime = 10, etime = 300;
  int time = etime - stime;
  int number, timestep, nodes=100*multi, posmin[time+1], edge = 10*multi, lnumber, snumber, allnumber, type;
  double xs, ys, zs, surfxs, ent, enmin, stress, volume,s;
  int group, count=0;
  double entropy[nodes], endif[nodes];
  int num[nodes];
  double dV = 40.78 * 40.78 * 40.78; // 1/A^3
  double dx = 40.78;

  fentro.open("./entropy.txt", ios::in);
  fliquid.open("./liquid.txt", ios::out);

  fliquid << scientific;

  while (getline(fentro,str)) {
      timestep = stoi(str);
      if (timestep > etime*1000) break;
      getline(fentro,str);
      streampos oldpos;
      for (int n=0; n < nodes; n++) {
        entropy[n] = 0;
        endif[n] = 0;
      }
      for (int n=0; n <= nodes; n++) {
        oldpos = fentro.tellg();
        fentro >> nodesid >> ent ;
        if (fentro.tellg() == -1) {
          fentro.clear();
          fentro.seekg(oldpos);
          getline(fentro,str);
          getline(fentro,str);
          break;
        } else {
          entropy[nodesid] = ent;
        }
      }
      for (int n=1+edge; n < nodes-edge; n++) {
        if (entropy[n]*entropy[n-1] != 0)
          endif[n] = entropy[n] - entropy[n-1];
      }
      int ii = timestep/1000 - stime;
      /*
      enmin = -0.005*dx;
      posmin[ii] = 0;
      for (int n=nodes-1; n >= 0; n--) {
        if (endif[n] < enmin) {
          posmin[ii] = n;
          break;
        }
      }
      */
      enmin = 0;
      if(ii >= 10) {
        start = posmin[ii-1] - 5;
      } else {
        start = 41;
      }
      if(start < 41) start = 41;
      if(start >= nodes) start = nodes;
      posmin[ii] = start;
      for (int n=start; n < start+10; n++) {
        if (endif[n] < enmin) {
          enmin = endif[n];
          posmin[ii] = n;
        }
      }
  }
  //return 0;
  fentro.close();

  for(int n=stime*1000; n<=etime*1000; n += 1000) {
    int ii = n/1000 - stime;

    input.open(pre + to_string(n) + end);

    for (int i=0; i < 9; i++) {
      getline(input,line);
      if (i == 1) {
        timestep = stoi(line);
        if (timestep != n) {
          cout << "Error in reading files\n";
          return 1;
        }
      }
      if (i == 3) number = stoi(line);
      if (ii == 0 && i == 3) allnumber = stoi(line);
    }
    cout << n/1000 << "\n";
    snumber = 0;
    for (int i=0; i < number; i++){
      input >> id >> type >> xs >> ys >> zs >> stress >> volume >> group >> s;

      if (xs*nodes > posmin[ii]) snumber++;
    }

    lnumber = allnumber - snumber;
    fliquid << n << " " << posmin[ii] << " " << lnumber <<"\n";
    input.close();
  }

  fliquid.close();

  return 0;

}
