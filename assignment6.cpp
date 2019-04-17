#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

using namespace std;
#define BigM -1
vector<vector<int> > comb_ind;

template <typename T>
void print2D(vector<vector<T> > mat) {
  int n_rows = mat.size();
  int n_cols = mat[0].size();
  for (int i = 0; i < n_rows; i++) {
    for (int j = 0; j < n_cols; j++) cout << mat[i][j] << "\t";
    cout << endl;
  }
}

template <typename T>
void print1D(vector<T> vec) {
  for (int i = 0; i < vec.size(); i++) {
    cout << vec[i] << " ";
  }
  cout << endl;
}

void print_simplex_table(vector<vector<double> > mat,
                         pair<vector<int>, vector<int> > format) {
  int n_rows = mat.size() - 1;
  int n_cols = mat[0].size();
  for (int i = 0; i < format.first.size(); i++) {
    if (format.first[i] <= format.first.size())
      cout << "-x" << format.first[i] << "\t";
    else
      cout << "-Z" << (format.first[i] - format.first.size()) << "\t";
  }
  cout << "1";
  cout << endl;
  for (int i = 0; i < n_rows; i++) {
    for (int j = 0; j < n_cols; j++) cout << mat[i][j] << "\t";

    // cout<<mat[i][mat[0].size()-1]<<"\t";
    if (format.second[i] <= format.first.size())
      cout << "x" << format.second[i] << "\t";
    else
      cout << "Z" << (format.second[i] - format.first.size()) << "\t";
    cout << endl;
  }

  for (int i = 0; i < n_cols; i++) cout << mat[mat.size() - 1][i] << "\t";
  cout << endl;
}

void simplex_solver(vector<vector<double> > full_eqns,
                    vector<double> full_obj) {
  // Initialization
  int rows = full_eqns.size();
  int cols = full_eqns[0].size();
  int var1 = 0,
      var2 = 0;  // var1 counts the artificial vars & var2 counts the slack vars
  for (int i = 0; i < rows; i++) {
    if (full_eqns[i][cols - 2] == -1) {
      var2++;
    } else
      var1++;
  }
  var1 += var2;
  int n_cols = cols + var2 - 1;  // modified new row len
  vector<vector<double> > eqns_final(rows + 1, vector<double>(n_cols, 0.0));
  int k = cols - 2;
  for (int i = 0; i < rows; i++) {
    bool ddone = true;
    for (int j = 0; j < n_cols; j++) {
      if (j < (cols - 2))
        eqns_final[i][j] = full_eqns[i][j];
      else if (j == n_cols - 1)  // get the constant in
        eqns_final[i][j] = full_eqns[i][cols - 1];
      else {
        if ((full_eqns[i][cols - 2] == -1) && (j == k) && ddone) {
          eqns_final[i][j] = -1.0;
          k++;
          ddone = false;
        }
      }
    }
  }
  cout << "FULL EQUATIONS" << endl;
  print2D(full_eqns);
  vector<double> M_index;  /// var1 denotes the length along the row --->
  for (int i = 0; i < var1; i++) {
    if (full_eqns[i][cols - 2] == -1)
      M_index.push_back(BigM);
    else
      M_index.push_back(0.0);
  }
  cout << "M_INDEX" << endl;
  print1D<double>(M_index);
  for (int i = 0; i < var2 - 1; i++) {
    full_obj.push_back(
        0.0);  // fill up the remaining slots of objective functions
  }
  cout << "Eqns" << endl;
  print1D<double>(full_obj);
  print2D<double>(eqns_final);
  for (int i = 0; i < n_cols; i++) {
    for (int j = 0; j < rows; j++) {
      // cout << i << "," << j << endl;
      eqns_final[rows][i] +=
          M_index[j] * eqns_final[j][i];  // getting the final last row values
    }
    // if (i < n_cols - 1) eqns_final[rows][i] += full_obj[i];
  }
  cout << "Initial Table" << endl;
  print2D<double>(eqns_final);
  // return;
  vector<int> sol_nb(n_cols - 1, 0);
  vector<int> sol_b(rows, 0);
  for (int i = 0; i < sol_nb.size(); i++) sol_nb[i] = i + 1;
  for (int i = 0; i < sol_b.size(); i++) sol_b[i] = i + var1 + rows;

  pair<vector<int>, vector<int> > indices;
  indices = make_pair(sol_nb, sol_b);
  print_simplex_table(eqns_final, indices);

  // book-keeping
  int count1 = 0;
  vector<vector<vector<double> > > table1_cache;
  vector<pair<vector<int>, vector<int> > > index1_cache;
  vector<vector<double> > thetas1_cache;
  vector<double> pivot1_cache;
  int count2 = 0;
  vector<vector<vector<double> > > table2_cache;
  vector<pair<vector<int>, vector<int> > > index2_cache;
  vector<vector<double> > thetas2_cache;
  vector<pair<int, int> > pivot2_cache;
  table1_cache.push_back(eqns_final);
  index1_cache.push_back(indices);
  bool unbounded = true;
  bool finale = true;  // the last case to handle a zero case
  // iters
  while (1) {
    vector<double> thetas(sol_b.size(), 0.0);
    int obj_row = eqns_final.size() - 1;
    cout << "Computing Thetas" << endl;
    int min_pos = 0;
    bool done = true;
    double min = 1e6;
    for (int i = 0; i < eqns_final[0].size() - 1; i++) {
      if ((min > eqns_final[obj_row][i]) && (eqns_final[obj_row][i] < 0)) {
        min_pos = i;
        min = eqns_final[obj_row][min_pos];
        done = false;
      }
    }
    /*
    if (done && finale)
    {
        for(int i = 0; i < full_obj.size(); i++)
        {
            if(eqns_final[obj_row][i] == 0.0)
            {
                done = false;
                finale = true;
            }
        }

    }*/
    unbounded = true;
    for (int i = 0; i < eqns_final.size(); i++) {
      if (eqns_final[i][min_pos] > 0)
        unbounded = false;  // when all the rqtios are -ve
    }

    if (done || unbounded) {
      break;
    }

    else
      count1++;

    for (int i = 0; i < thetas.size(); i++)
      thetas[i] =
          eqns_final[i][eqns_final[0].size() - 1] / eqns_final[i][min_pos];
    cout << "MIN_POSE = " << min_pos << endl;
    cout << "THETAS: ";
    print1D(thetas);
    thetas1_cache.push_back(thetas);
    int min_theta_pose = 0;
    min = 1e6;
    for (int i = 0; i < thetas.size(); i++) {
      if (thetas[i] < 0) continue;
      if (min >= thetas[i]) {
        min_theta_pose = i;
        min = thetas[i];
      }
    }

    cout << "MIN_theta_POSE = " << min_theta_pose << endl;
    cout << "Pivot Element = " << eqns_final[min_theta_pose][min_pos] << endl;
    double pivot = eqns_final[min_theta_pose][min_pos];
    pivot1_cache.push_back(pivot);
    pair<int, int> pv(min_theta_pose, min_pos);
    int temp = 0;
    temp = indices.first[pv.second];
    indices.first[pv.second] = indices.second[pv.first];
    indices.second[pv.first] = temp;
    index1_cache.push_back(indices);
    for (int i = 0; i < eqns_final.size(); i++) {
      for (int j = 0; j < eqns_final[0].size(); j++) {
        if (i == pv.first || j == pv.second) continue;
        double p = 0.0, q = 0.0, r = 0.0, s = 0.0;
        p = pivot;
        s = eqns_final[i][j];
        r = eqns_final[pv.first][j];
        q = eqns_final[i][pv.second];
        eqns_final[i][j] = (p * s - q * r) / p;
      }
    }
    for (int i = 0; i < eqns_final[0].size(); i++) {
      if (i == pv.second) continue;
      eqns_final[pv.first][i] /= pivot;
    }
    for (int i = 0; i < eqns_final.size(); i++) {
      if (i == pv.first) continue;
      eqns_final[i][pv.second] /= -pivot;
    }
    eqns_final[pv.first][pv.second] = 1 / pivot;
    table1_cache.push_back(eqns_final);
    print_simplex_table(eqns_final, indices);
  }
  cout << "DONE FINAL" << endl;
  print_simplex_table(eqns_final, indices);
  print1D(indices.first);
  print1D(indices.second);
  print1D(full_obj);
  // cout<<"MUHAHAHAHA"<<endl;
  int num_var = indices.first.size();
  for (int i = 0; i < eqns_final[0].size() - 1; i++) {
    if (indices.first[i] > num_var) {
      for (int j = 0; j < eqns_final.size(); j++) {
        // cout<<i<<","<<j<<":"<<eqns_final.size()<<endl;
        eqns_final[j][i] = 0.0;
      }
    }
  }
  for (int i = 0; i < eqns_final[0].size(); i++) {
    if (indices.first[i] <= num_var || i == eqns_final[0].size() - 1) {
      double temp = 0.0;
      for (int j = 0; j < eqns_final.size() - 1; j++) {
        // cout<<i<<","<<j<<":"<<eqns_final.size()<<endl;
        if (indices.second[j] <= indices.second.size()) {
          temp += eqns_final[j][i] * -(full_obj[indices.second[j] - 1]);
          // cout<<"TEMP"<<temp<<endl;
          // cout<<"OBJ"<<full_obj[indices.second[j]-1]<<endl;
          // cout<<"INDEX"<<indices.second[j]<<endl;
        }
      }
      double k = 0.0;
      if (indices.first[i] <= indices.second.size())
        k = full_obj[indices.first[i] - 1];
      eqns_final[eqns_final.size() - 1][i] = temp + k;
    }
  }
  cout << "PHASE 2 START" << endl;
  print_simplex_table(eqns_final, indices);
  full_eqns = eqns_final;
  print_simplex_table(full_eqns, indices);
  int n = full_eqns.size();
  int r = full_eqns[0].size() - 1;
  table2_cache.push_back(full_eqns);
  index2_cache.push_back(indices);
  while (1) {
    vector<double> thetas(n, 0.0);
    int obj_row = full_eqns.size() - 1;
    cout << "Computing Thetas" << endl;
    int min_pos = 0;
    double min = 1e6;
    bool done = true;
    for (int i = 0; i < full_eqns[0].size() - 1; i++) {
      if ((min > full_eqns[obj_row][i]) && full_eqns[obj_row][i] < 0) {
        min_pos = i;
        min = full_eqns[obj_row][min_pos];
        // cout<<"YAAY"<<endl;
        done = false;
      }
    }
    // cout<<"MIN MIN"<<min_pos;
    unbounded = true;
    for (int i = 0; i < full_eqns.size(); i++) {
      if (full_eqns[i][min_pos] > 0) unbounded = false;
    }

    if (done || unbounded) {
      break;
    }

    else
      count2++;

    for (int i = 0; i < thetas.size(); i++)
      thetas[i] = full_eqns[i][r] / full_eqns[i][min_pos];
    cout << "MIN_POSE = " << min_pos << endl;
    cout << "THETAS: ";
    print1D(thetas);
    thetas2_cache.push_back(thetas);
    int min_theta_pose = 0;
    double min_theta = 1e6;
    for (int i = 0; i < thetas.size(); i++) {
      if (thetas[i] < 0) continue;
      if (min_theta > thetas[i]) {
        // cout<<"YAAY"<<endl;
        cout<<"MIN = "<<min<<endl;
        min_theta_pose = i;
        min_theta = thetas[min_theta_pose];
      }
    }

    cout << "MIN_theta_POSE = " << min_theta_pose << endl;
    cout << "Pivot Element = " << full_eqns[min_theta_pose][min_pos] << endl;
    double pivot = full_eqns[min_theta_pose][min_pos];
    pair<int, int> pv(min_theta_pose, min_pos);
    pivot2_cache.push_back(pv);

    int temp = 0;
    temp = indices.first[pv.second];
    indices.first[pv.second] = indices.second[pv.first];
    indices.second[pv.first] = temp;
    index2_cache.push_back(indices);
    for (int i = 0; i < full_eqns.size(); i++) {
      for (int j = 0; j < full_eqns[0].size(); j++) {
        if (i == pv.first || j == pv.second) continue;
        double p = 0.0, q = 0.0, r = 0.0, s = 0.0;
        p = pivot;
        s = full_eqns[i][j];
        r = full_eqns[pv.first][j];
        q = full_eqns[i][pv.second];
        full_eqns[i][j] = (p * s - q * r) / p;
      }
    }
    for (int i = 0; i < full_eqns[0].size(); i++) {
      if (i == pv.second) continue;
      full_eqns[pv.first][i] /= pivot;
    }
    for (int i = 0; i < full_eqns.size(); i++) {
      if (i == pv.first) continue;
      full_eqns[i][pv.second] /= -pivot;
    }
    full_eqns[pv.first][pv.second] = 1 / pivot;
    table2_cache.push_back(full_eqns);
    print_simplex_table(full_eqns, indices);
    // return;
  }

  while (1) {
    cout << "a: Phase1 iters \nb: Phase1 BFS\n"
            "c: Phase1 objective value \nd: Phase1 Optimal Solution \n"
            "e: Phase2 iters \nf: Phase2 BFS\n"
            "g: Phase2 objective value \nh: Phase2 Optimal Solution\n"
            "i: most negative zj-cj\n"
            "x: To end"
         << endl;
    char ch;
    cin >> ch;
    if (ch == 'x') break;
    switch (ch) {
      case 'a': {
        cout << count1 << endl;
        break;
      }

      case 'b': {
        cout << "Enter Iter Number <=  " << count1 << endl;
        int val;
        cin >> val;
        for (int j = 0; j < index1_cache[val].first.size(); j++) {
          if (index1_cache[val].first[j] <= index1_cache[val].first.size())
            cout << "x" << index1_cache[val].first[j] << ":\t";
          else
            cout << "Z" << (index1_cache[val].first[j] -
                            index1_cache[val].first.size())
                 << ":\t";
          cout << "0.0"
               << "\t";
        }
        cout << endl;

        for (int j = 0; j < index1_cache[val].second.size(); j++) {
          if (index1_cache[val].second[j] <= index1_cache[val].first.size())
            cout << "x" << index1_cache[val].second[j] << ":\t";
          else
            cout << "Z" << (index1_cache[val].second[j] -
                            index1_cache[val].first.size())
                 << ":\t";
          cout << table1_cache[val][j][eqns_final[0].size() - 1] << "\t";
        }
        cout << endl;

        break;
      }

      case 'c': {
        cout << "Enter ith iter" << endl;
        int val;
        cin >> val;
        cout << "Optimal Value: " << table1_cache[val][eqns_final.size() - 1]
                                                 [eqns_final[0].size() - 1]
             << endl;
        break;
      }

      case 'd': {
        cout << "Optimal Value: " << table1_cache[table1_cache.size() - 1]
                                                 [eqns_final.size() - 1]
                                                 [eqns_final[0].size() - 1]
             << endl;
        for (int j = 0; j < index1_cache[index1_cache.size() - 1].first.size();
             j++) {
          if (index1_cache[index1_cache.size() - 1].first[j] <=
              index1_cache[index1_cache.size() - 1].first.size())
            cout << "x" << index1_cache[index1_cache.size() - 1].first[j]
                 << ":\t";
          else
            cout << "Z" << (index1_cache[index1_cache.size() - 1].first[j] -
                            index1_cache[index1_cache.size() - 1].first.size())
                 << ":\t";
          cout << "0.0"
               << "\t";
        }
        cout << endl;

        for (int j = 0; j < index1_cache[index1_cache.size() - 1].second.size();
             j++) {
          if (index1_cache[index1_cache.size() - 1].second[j] <=
              index1_cache[index1_cache.size() - 1].first.size())
            cout << "x" << index1_cache[index1_cache.size() - 1].second[j]
                 << ":\t";
          else
            cout << "Z" << (index1_cache[index1_cache.size() - 1].second[j] -
                            index1_cache[index1_cache.size() - 1].first.size())
                 << ":\t";
          cout << table1_cache[index1_cache.size() - 1][j]
                              [eqns_final[0].size() - 1]
               << "\t";
        }
        cout << endl;

        break;
      }

      case 'e': {
        cout << count2 << endl;
        break;
      }

      case 'f': {
        cout << "Enter Iter Number <=  " << count2 << endl;
        int val;
        cin >> val;
        for (int j = 0; j < index2_cache[val].first.size(); j++) {
          if (index2_cache[val].first[j] <= index2_cache[val].first.size())
            cout << "x" << index2_cache[val].first[j] << ":\t";
          else
            cout << "Z" << (index2_cache[val].first[j] -
                            index2_cache[val].first.size())
                 << ":\t";
          cout << "0.0"
               << "\t";
        }
        cout << endl;

        for (int j = 0; j < index2_cache[val].second.size(); j++) {
          if (index2_cache[val].second[j] <= index2_cache[val].first.size())
            cout << "x" << index2_cache[val].second[j] << ":\t";
          else
            cout << "Z" << (index2_cache[val].second[j] -
                            index2_cache[val].first.size())
                 << ":\t";
          cout << table2_cache[val][j][eqns_final[0].size() - 1] << "\t";
        }
        cout << endl;

        break;
      }

      case 'g': {
        cout << "Enter ith iter" << endl;
        int val;
        cin >> val;
        cout << "Optimal Value: " << table2_cache[val][eqns_final.size() - 1]
                                                 [eqns_final[0].size() - 1]
             << endl;
        break;
      }

      case 'h': {
        cout << "Optimal Value: " << table2_cache[table2_cache.size() - 1]
                                                 [eqns_final.size() - 1]
                                                 [eqns_final[0].size() - 1]
             << endl;
        for (int j = 0; j < index2_cache[index2_cache.size() - 1].first.size();
             j++) {
          if (index2_cache[index2_cache.size() - 1].first[j] <=
              index2_cache[index2_cache.size() - 1].first.size())
            cout << "x" << index2_cache[index2_cache.size() - 1].first[j]
                 << ":\t";
          else
            cout << "Z" << (index2_cache[index2_cache.size() - 1].first[j] -
                            index2_cache[index2_cache.size() - 1].first.size())
                 << ":\t";
          cout << "0.0"
               << "\t";
        }
        cout << endl;

        for (int j = 0; j < index2_cache[index2_cache.size() - 1].second.size();
             j++) {
          if (index2_cache[index2_cache.size() - 1].second[j] <=
              index2_cache[index2_cache.size() - 1].first.size())
            cout << "x" << index2_cache[index2_cache.size() - 1].second[j]
                 << ":\t";
          else
            cout << "Z" << (index2_cache[index2_cache.size() - 1].second[j] -
                            index2_cache[index2_cache.size() - 1].first.size())
                 << ":\t";
          cout << table2_cache[index2_cache.size() - 1][j]
                              [eqns_final[0].size() - 1]
               << "\t";
        }
        cout << endl;
        break;
      }

      case 'i': {
        cout << "Enter Iter Number <=  " << count2 << endl;
        int val;
        cin >> val;

        cout << table2_cache[val][table2_cache[0].size() - 1]
                            [pivot2_cache[val].second];
        break;
      }

      default: { cout << "Invalid Option"; }
    }
  }
}

int main(int argc, char const *argv[]) {
  cout << "Enter the number of Equations" << endl;
  int n, r;
  cin >> n;
  cout << "Enter the number of variables" << endl;
  cin >> r;
  cout << "a1x1 + a2x2 + ... + arxr (>/</=)(specify when asked -1/+1/0) C"
       << endl;
  vector<vector<double> > eqns;
  // vector<vector<double> > eqns_only(n, vector<double>(r + 1, 0.0));
  // vector<vector<double> > eqns_format(n, vector<double>((n + r + 1), 0.0));
  for (int i = 0; i < n; i++) {
    vector<double> arr(2 + r, 0);
    for (int j = 0; j < r; j++) {
      cin >> arr[j];
      // eqns_format[i][j] = arr[j];
      // eqns_only[i][j] = arr[j];
    }
    cout << "enter sign" << endl;
    cin >> arr[r];
    // eqns_format[i][i + r] = arr[r];
    cout << "enter C" << endl;
    cin >> arr[r + 1];
    // eqns_format[i][n + r] = arr[r + 1];
    // eqns_only[i][r] = arr[r + 1];
    cout << "Got Eqn: " << i << endl;
    eqns.push_back(arr);
  }
  cout << "Equations Data" << endl;
  print2D<double>(eqns);

  cout << "Enter Final Function coefficients for maximize(1) or minimize(-1)"
       << endl;
  int eq_flag;
  cin >> eq_flag;
  cout << "Enter Coefficients" << endl;
  vector<double> optimal_eq;
  vector<double> opt_full_splx(n + r, 0.0);
  for (int i = 0; i < r; i++) {
    double temp = 0.0;
    cin >> temp;
    opt_full_splx[i] = -temp;
    optimal_eq.push_back(-temp);
  }
  optimal_eq.push_back(0.0);
  print1D(opt_full_splx);
  simplex_solver(eqns, optimal_eq);
  return 0;
}
