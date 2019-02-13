#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

using namespace std;

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
  int r = format.first.size();
  for (int i = 0; i < r; i++) {
    if (format.first[i] <= r)
      cout << "-x" << format.first[i] << "\t";
    else
      cout << "-Z" << (format.first[i] - r) << "\t";
  }
  cout << "1";
  cout << endl;
  // top header printing ^
  for (int i = 0; i < n_rows; i++) {
    for (int j = 0; j < n_cols; j++) cout << mat[i][j] << "\t";

    // cout<<mat[i][mat[0].size()-1]<<"\t";
    if (format.second[i] <= r)
      cout << "x" << format.second[i] << "\t";
    else
      cout << "Z" << (format.second[i] - r) << "\t";
    cout << endl;
  }

  for (int i = 0; i < n_cols; i++)
    cout << mat[mat.size() - 1][i] << "\t";  // for the constants
  cout << endl;
}

void simplex_solver(vector<vector<double> > full_eqns,
                    vector<double> full_obj) {
  // Initialization
  int n = full_eqns.size();         // rows
  int r = full_eqns[0].size() - 1;  // cols
  pair<vector<int>, vector<int> > indices;
  vector<int> sol_nb(r, 0);
  vector<int> sol_b(n, 0);
  for (int i = 0; i < sol_nb.size(); i++) sol_nb[i] = i + 1;  // the zero vars
  for (int i = 0; i < sol_b.size(); i++)
    sol_b[i] = i + 1 + r;  // the non zero vars
  indices = make_pair(sol_nb, sol_b);
  full_eqns.push_back(full_obj);
  cout << "Initial Table" << endl;
  print_simplex_table(full_eqns, indices);
  // book-keeping
  int count = 0;
  vector<vector<vector<double> > > table_cache;
  vector<pair<vector<int>, vector<int> > > index_cache;
  vector<vector<double> > thetas_cache;
  table_cache.push_back(full_eqns);
  index_cache.push_back(indices);
  bool unbounded = true;
  // iters
  while (1) {
    vector<double> thetas(n, 0.0);
    int obj_row = full_eqns.size() - 1;
    cout << "Computing Thetas" << endl;
    int min_pos = 0;
    bool done = true;
    for (int i = 0; i < full_obj.size() - 1; i++) {
      if (full_eqns[obj_row][min_pos] >= full_eqns[obj_row][i] &&
          full_eqns[obj_row][i] < 0) {
        min_pos = i;
        done = false;
      }
    }
    unbounded = true;
    for (int i = 0; i < full_eqns.size(); i++) {
      if (full_eqns[i][min_pos] > 0) unbounded = false;
    }

    if (done || unbounded) {
      break;
    }

    else
      count++;

    for (int i = 0; i < thetas.size(); i++)
      thetas[i] = full_eqns[i][r] / full_eqns[i][min_pos];
    cout << "MIN_POSE = " << min_pos << endl;
    cout << "THETAS: ";
    print1D(thetas);
    thetas_cache.push_back(thetas);
    int min_theta_pose = 0;

    for (int i = 0; i < thetas.size(); i++) {
      if (thetas[i] < 0) continue;
      if (thetas[min_theta_pose] >= thetas[i]) {
        min_theta_pose = i;
      }
    }

    cout << "MIN_theta_POSE = " << min_theta_pose << endl;
    cout << "Pivot Element = " << full_eqns[min_theta_pose][min_pos] << endl;
    double pivot = full_eqns[min_theta_pose][min_pos];
    pair<int, int> pv(min_theta_pose, min_pos);

    int temp = 0;
    temp = indices.first[pv.second];
    indices.first[pv.second] = indices.second[pv.first];
    indices.second[pv.first] = temp;
    index_cache.push_back(indices);
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
    table_cache.push_back(full_eqns);
    print_simplex_table(full_eqns, indices);
  }
  cout << "DONE" << endl;
  cout << "a: Number of Iters \nb: Non Basic Vals \nc: Basic with Min_ratio\n"
          "d: Simplex Table \ne: Optimal Solution \nx: To end"
       << endl;
  while (1) {
    char ch;
    cin >> ch;
    if (ch == 'x') break;
    switch (ch) {
      case 'a': {
        cout << "Iterations: " << count << endl;
        break;
      }

      case 'b': {
        cout << "Enter ith iter" << endl;
        int i;
        cin >> i;
        i = i + 1;
        cout << "Non Basic Solutions" << endl;
        // print1D<int>(index_cache[i-1].first);
        for (int j = 0; j < index_cache[i - 1].first.size(); j++) {
          if (index_cache[i - 1].first[j] <= index_cache[i - 1].first.size())
            cout << "-x" << index_cache[i - 1].first[j] << "\t";
          else
            cout << "-Z" << (index_cache[i - 1].first[j] -
                             index_cache[i - 1].first.size())
                 << "\t";
        }
        cout << endl;
        cout << "COMPUTE :" << endl;
        for (int j = 0; j < full_eqns.size(); j++) {
          cout << table_cache[i - 1][full_eqns[0].size() - 1][j] << "\t";
        }
        cout << endl;
        break;
      }

      case 'c': {
        cout << "Enter ith iter" << endl;
        int i;
        cin >> i;
        i = i + 1;
        cout << "Basic Solutions" << endl;
        print1D<int>(index_cache[i - 1].second);
        for (int j = 0; j < index_cache[i - 1].second.size(); j++) {
          if (index_cache[i - 1].second[j] <= index_cache[i - 1].first.size())
            cout << "x" << index_cache[i - 1].second[j] << "\t";
          else
            cout << "Z" << (index_cache[i - 1].second[j] -
                            index_cache[i - 1].first.size())
                 << "\t";
        }
        cout << endl;
        cout << "Thetas" << endl;
        print1D<double>(thetas_cache[i - 2]);
        cout << "MIN RATIO ";
        int min_theta_pose = 0;
        for (int i = 0; i < thetas_cache[i - 2].size(); i++) {
          if (thetas_cache[i - 2][i] < 0) continue;
          if (thetas_cache[i - 2][min_theta_pose] >= thetas_cache[i - 2][i]) {
            min_theta_pose = i;
          }
        }
        cout << thetas_cache[i - 2][min_theta_pose] << endl;
        break;
      }

      case 'd': {
        cout << "Enter i for ith table" << endl;
        int i;
        cin >> i;
        print_simplex_table(table_cache[i - 1], index_cache[i - 1]);
        cout << "Basic Feasible Solutions" << endl;
        print1D<int>(index_cache[i - 1].second);
        cout << "Thetas" << endl;
        print1D<double>(thetas_cache[i - 2]);
        break;
      }

      case 'e': {
        bool infeasible = false;
        for (int j = 0; j < indices.second.size(); j++) {
          if (full_eqns[j][full_eqns[0].size() - 1] < 0) infeasible = true;
        }
        if (infeasible) {
          cout << "Solution is infeasible" << endl;
        }

        if (unbounded) {
          cout << "Solution is Unbounded" << endl;
        }

        if (!infeasible && !unbounded) {
          cout << "Optimal Value: " << full_eqns[full_eqns.size() - 1]
                                                [full_eqns[0].size() - 1]
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
            cout << full_eqns[j][full_eqns[0].size() - 1] << "\t";
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
  vector<vector<double> > eqns_only(n, vector<double>(r + 1, 0.0));
  // vector<vector<double> > eqns_format(n, vector<double>((n + r + 1), 0.0));
  for (int i = 0; i < n; i++) {
    vector<double> arr(2 + r, 0);
    for (int j = 0; j < r; j++) {
      cin >> arr[j];
      // eqns_format[i][j] = arr[j];
      eqns_only[i][j] = arr[j];
    }
    cout << "enter sign" << endl;
    cin >> arr[r];
    // eqns_format[i][i + r] = arr[r];
    cout << "enter C" << endl;
    cin >> arr[r + 1];
    // eqns_format[i][n + r] = arr[r + 1];
    eqns_only[i][r] = arr[r + 1];
    cout << "Got Eqn: " << i << endl;
    eqns.push_back(arr);
  }
  cout << "Equations Data" << endl;
  print2D<double>(eqns_only);
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
  simplex_solver(eqns_only, optimal_eq);
  return 0;
}
