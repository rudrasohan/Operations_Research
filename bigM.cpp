#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

using namespace std;
#define BigM 1e4
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
  int var1 = 0, var2 = 0;//var1 counts the artificial vars & var2 counts the slack vars
  for (int i = 0; i < rows; i++) {
    if (full_eqns[i][cols - 2] == -1) {
      var2++;
    } else
      var1++;
  }
  var1 += var2; 
  int n_cols = cols + var2 - 1; //modified new row len
  vector<vector<double> > eqns_final(rows + 1, vector<double>(n_cols, 0.0));
  int k = cols - 2;
  for (int i = 0; i < rows; i++) {
    bool ddone = true;
    for (int j = 0; j < n_cols; j++) {
      if (j < (cols - 2))
        eqns_final[i][j] = full_eqns[i][j];
      else if (j == n_cols - 1) //get the constant in
        eqns_final[i][j] = full_eqns[i][cols - 1];
      else {
        if ((full_eqns[i][cols - 2] == -1) &&
            (j == k) && ddone) {
          eqns_final[i][j] = -1.0;
          k++;
          ddone = false;
        }
      }
    }
  }
  vector<double> M_index; ///var1 denotes the length along the row --->
  for (int i = 0; i < var1; i++) {
    M_index.push_back(-BigM);
  }
  cout << "M_INDEX" << endl;
  print1D<double>(M_index);
  for (int i = 0; i < var2 - 1; i++) {
    full_obj.push_back(0.0); //fill up the remaining slots of objective functions
  }
  cout << "Eqns" << endl;
  print1D<double>(full_obj);
  print2D<double>(eqns_final);
  for (int i = 0; i < n_cols; i++) {
    for (int j = 0; j < rows; j++) {
      //cout << i << "," << j << endl;
      eqns_final[rows][i] += M_index[j] * eqns_final[j][i]; //getting the final last row values
    }
    if (i < n_cols - 1) eqns_final[rows][i] += full_obj[i];
  }
  cout << "Initial Table" << endl;
  print2D<double>(eqns_final);
  vector<int> sol_nb(n_cols - 1, 0);
  vector<int> sol_b(rows, 0);
  for (int i = 0; i < sol_nb.size(); i++) sol_nb[i] = i + 1;
  for (int i = 0; i < sol_b.size(); i++) sol_b[i] = i + var1 + rows;

  pair<vector<int>, vector<int> > indices;
  indices = make_pair(sol_nb, sol_b);
  print_simplex_table(eqns_final, indices);

  // book-keeping
  int count = 0;
  vector<vector<vector<double> > > table_cache;
  vector<pair<vector<int>, vector<int> > > index_cache;
  vector<vector<double> > thetas_cache;
  vector<double> pivot_cache;
  table_cache.push_back(eqns_final);
  index_cache.push_back(indices);
  bool unbounded = true;
  bool finale = true; //the last case to handle a zero case
  // iters
  while (1) {
    vector<double> thetas(sol_b.size(), 0.0);
    int obj_row = eqns_final.size() - 1;
    cout << "Computing Thetas" << endl;
    int min_pos = 0;
    bool done = true;
    double min = 1e6;
    for (int i = 0; i < full_obj.size(); i++) {
      if ((min >= eqns_final[obj_row][i]) && (eqns_final[obj_row][i] < 0)) {
        min_pos = i;
        min = eqns_final[obj_row][min_pos];
        done = false;
      }
    }
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
        
    }
    unbounded = true;
    for (int i = 0; i < eqns_final.size(); i++) {
      if (eqns_final[i][min_pos] > 0) unbounded = false; //when all the rqtios are -ve
    }

    if (done || unbounded) {
      break;
    }

    else
      count++;

    for (int i = 0; i < thetas.size(); i++)
      thetas[i] =
          eqns_final[i][eqns_final[0].size() - 1] / eqns_final[i][min_pos];
    cout << "MIN_POSE = " << min_pos << endl;
    cout << "THETAS: ";
    print1D(thetas);
    thetas_cache.push_back(thetas);
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
    pivot_cache.push_back(pivot);
    pair<int, int> pv(min_theta_pose, min_pos);
    int temp = 0;
    temp = indices.first[pv.second];
    indices.first[pv.second] = indices.second[pv.first];
    indices.second[pv.first] = temp;
    index_cache.push_back(indices);
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
      if (i == pv.second)
        continue;
      eqns_final[pv.first][i] /= pivot;
    }
    for (int i = 0; i < eqns_final.size(); i++) {
      if (i == pv.first)
        continue;
      eqns_final[i][pv.second] /= -pivot;
    }
    eqns_final[pv.first][pv.second] = 1/pivot; 
    table_cache.push_back(eqns_final);
    print_simplex_table(eqns_final, indices);
  }
  cout << "DONE" << endl;

  while (1) {
    cout << "a: Non Basic Vals \nb: Min_ratio\n"
            "c: Simplex Table \nd: Optimal Solution \nx: To end"
         << endl;
    char ch;
    cin >> ch;
    if (ch == 'x') break;
    switch (ch) {
      case 'a': {
        cout << "Enter ith iter" << endl;
        int i, k;
        cin >> i;
        cout << "Enter jth column" << endl;
        cin >> k;

        for (int j = 0; j < eqns_final.size(); j++) {
          cout << table_cache[i][j][k] << "\t";
        }
        cout << endl;
        break;
      }

      case 'b': {
        cout << "Enter ith iter" << endl;
        int i;
        cin >> i;
        cout << "Ratios:" << endl;
        print1D<double>(thetas_cache[i - 1]);
        cout << "PIVOT:" << pivot_cache[i - 1] << endl;
        break;
      }

      case 'c': {
        cout << "Enter ith iter" << endl;
        int i;
        cin >> i;
        print_simplex_table(table_cache[i], index_cache[i]);
        break;
      }

      case 'd': {
        bool infeasible = false;
        for (int j = 0; j < indices.second.size(); j++) {
          if (eqns_final[j][eqns_final[0].size() - 1] < 0) infeasible = true;
        }
        if (infeasible) {
          cout << "Solution is infeasible" << endl;
        }

        if (unbounded) {
          cout << "Solution is Unbounded" << endl;
        }

        if (!infeasible && !unbounded) {
          cout << "Optimal Value: " << eqns_final[eqns_final.size() - 1]
                                                 [eqns_final[0].size() - 1]
               << endl;
          for (int j = 0; j < indices.first.size(); j++) {
            if (indices.first[j] <= indices.first.size())
              cout << "x" << indices.first[j] << ":\t";
            else
              cout << "Z" << (indices.first[j] - indices.first.size()) << ":\t";
            cout << "0.0"
                 << "\t";
          }
          cout << endl;

          for (int j = 0; j < indices.second.size(); j++) {
            if (indices.second[j] <= indices.first.size())
              cout << "x" << indices.second[j] << ":\t";
            else
              cout << "Z" << (indices.second[j] - indices.first.size())
                   << ":\t";
            cout << eqns_final[j][eqns_final[0].size() - 1] << "\t";
          }
          cout << endl;
        }
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
  //vector<vector<double> > eqns_only(n, vector<double>(r + 1, 0.0));
  //vector<vector<double> > eqns_format(n, vector<double>((n + r + 1), 0.0));
  for (int i = 0; i < n; i++) {
    vector<double> arr(2 + r, 0);
    for (int j = 0; j < r; j++) {
      cin >> arr[j];
      //eqns_format[i][j] = arr[j];
      //eqns_only[i][j] = arr[j];
    }
    cout << "enter sign" << endl;
    cin >> arr[r];
    //eqns_format[i][i + r] = arr[r];
    cout << "enter C" << endl;
    cin >> arr[r + 1];
    //eqns_format[i][n + r] = arr[r + 1];
    //eqns_only[i][r] = arr[r + 1];
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
