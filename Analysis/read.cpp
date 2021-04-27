// g++ -std=c++11 ./read.cpp -o ./read;
#include <iostream>
#include <fstream>
using namespace std;

int main (int argc, char* argv[]) {
  string line, pre="./Au_p0.1_i17500_", end=".dump";
  ofstream output, fpress, fgroup, fentro, fdense, fvolume, fsurf;
  ifstream input;
  int id;
  int multi = 1;
  int number, timestep, nodes=100*multi, type, stime = 10, etime = 300;
  double xs, ys, zs, surfxs;
  int group, surN=30*30*2, count=0;
  double surX[surN];
  double press, stress, volume, s;
  double pressure[nodes], entropy[nodes], nu[nodes];
  int num[nodes];
  double dV = 122.3 * 122.3 * 40.78; // 1/A^3

  // output.open("./output.txt", ios::out);
  fpress.open("./pressure.txt", ios::out);
  fentro.open("./entropy.txt", ios::out);
  fdense.open("./density.txt", ios::out);
  fvolume.open("./volume.txt", ios::out);
  fsurf.open("./surface.txt", ios::out);
  //fgroup.open("./group.txt", ios::out);
  // output << scientific;
  fpress << scientific;
  fentro << scientific;
  fdense << scientific;
  fvolume << scientific;
  fsurf << scientific;
  //fgroup << scientific;

  for(int n=stime*1000; n<=etime*1000; n += 1000) {

    input.open(pre + to_string(n) + end);
    for (int i=0; i<nodes; i++) {
      num[i] = 0;
      pressure[i] = 0.0;
      entropy[i] = 0.0;
      nu[i] = 0.0;
    }
    for (int i=0; i < surN; i++) {
      surX[i] = 1;
    }

    for (int i=0; i < 9; i++) {
      getline(input,line);
      if (i == 1) {
        timestep = stoi(line);
        if (timestep != n) {
          cout << "Error in reading files\n";
          return 1;
        }
        fpress << timestep << "\n--------------------------------------------------------------------\n";
        fentro << timestep << "\n--------------------------------------------------------------------\n";
        fdense << timestep << "\n--------------------------------------------------------------------\n";
        fvolume << timestep << "\n--------------------------------------------------------------------\n";
        // fsurf << timestep << "\n--------------------------------------------------------------------\n";
        //fgroup << timestep << "\n--------------------------------------------------------------------\n";
      }
      if (i == 3) number = stoi(line);
    }
    cout << n/1000 << "\n";
    for (int i=0; i < number; i++){
      input >> id >> type >> xs >> ys >> zs >> stress >> volume >> group >> s;
      // output << id << " " << xs << " " << stress/volume << " " << group << "\n";
      int gg = 1;
      // if (n>=8000) gg = 2;
      if (group == gg) {
         //cout << count << "\n";
        if (count == 0) {
          surX[0] = xs;
          count ++;
        } else if (count == 1) {
          if (xs < surX[0]) {
            surX[1] = surX[0];
            surX[0] = xs;
          } else {
            surX[1] = xs;
          }
          count++;
        } else if (count == surN) {
          int myp = 0;
          if (xs < surX[surN-1]) {
            for (int k=0; k < surN-1; k++) {
              if (surX[k] <= xs && surX[k+1] > xs) {
                myp = k + 1;
                break;
              }
            }
            for (int k=surN-1; k > myp; k--) {
              surX[k] = surX[k-1];
            }
            surX[myp] = xs;
          }
        } else {
          int myp = 0;
          if (xs > surX[count-1]) {
            myp = count;
          }
          else {
            for (int k=0; k < count - 1; k++) {
              if (surX[k] <= xs && surX[k+1] > xs) {
                myp = k + 1;
                break;
              }
            }
          }
          for (int k=count; k > myp; k--) {
              surX[k] = surX[k-1];
            }
          surX[myp] = xs;
          count++;
        }
      }
      for (int j=0; j < nodes; j++) {
        if(xs >= double(j)/nodes && xs < double(j+1)/nodes) {
          num[j] += 1;
          pressure[j] += stress/volume;
          entropy[j] += s;
          nu[j] += volume;
          break;
        }
      }
    }
    surfxs = 0;
    for (int i=0; i < surN; i++) {
      // cout << i << " " << surX[i] << "\n";
      surfxs += surX[i];
    }
    surfxs /= surN;
    fsurf << n/1000 << " " << surfxs << "\n";
    for (int i=0; i < nodes; i++) {
      if (num[i] != 0 && num[i] >= 0.1 * 4000 * 3*3) {
        fpress << i << " " << pressure[i]/num[i] << "\n";
        fentro << i << " " << entropy[i]/num[i] << "\n";
        fdense << i << " " << num[i] << "\n";
        fvolume << i << " " << nu[i]/num[i] << "\n";
      }
    }

    fpress << "\n";
    fentro << "\n";
    fdense << "\n";
    fvolume << "\n";
    // fsurf << "\n";
    //fgroup << "\n";
    input.close();
  }
  // cout << id << " " << xs << " " << ys << " " << zs << " " << stress << " " << volume << " " << group;
  // output << line;
  // cout << scientific;
  // cout << 10000.0;
  // output.close();
  fpress.close();
  fentro.close();
  fdense.close();
  fvolume.close();
  fsurf.close();
  //fgroup.close();
  return 0;
}
