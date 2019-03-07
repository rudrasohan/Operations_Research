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

double simplex_solver(vector<vector<double> > full_eqns,
                      vector<double> full_obj, int status) {
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
  vector<pair<int, int> > pivot_cache;
  table_cache.push_back(full_eqns);
  index_cache.push_back(indices);
  bool unbounded = true;
  // iters
  while (1) {
    vector<double> thetas(full_eqns[0].size() - 1, 0.0);
    vector<double> thetas_no_mod(full_eqns[0].size() - 1, 0.0);
    int obj_row = full_eqns.size() - 1;
    int last_col = full_eqns[0].size() - 1;
    cout << "Computing Thetas" << endl;
    int min_pos = 0;
    bool done = true;
    /*
    for (int i = 0; i < full_obj.size() - 1; i++) {
      if (full_eqns[obj_row][min_pos] >= full_eqns[obj_row][i] &&
          full_eqns[obj_row][i] < 0) {
        min_pos = i;
        done = false;
      }
    }*/
    for (int i = 0; i < full_eqns.size() - 1; i++) {
      if (full_eqns[min_pos][last_col] >= full_eqns[i][last_col] &&
          full_eqns[i][last_col] < 0) {
        min_pos = i;
        done = false;
      }
    } /*
     unbounded = true;
     for (int i = 0; i < full_eqns[0].size(); i++) {
       if (full_eqns[i][min_pos] > 0) unbounded = false;
     }*/
    // cout<<"MIN ELE COL LOC =  "<<min_pos<<endl;
    // cout<<"MIN ELE COL"<<full_eqns[min_pos][last_col]<<endl;
    unbounded = true;
    for (int i = 0; i < full_eqns[0].size(); i++) {
      if (full_eqns[min_pos][i] < 0) unbounded = false;
    }

    if (done || unbounded) {
      if (done) unbounded = false;
      break;
    }

    else
      count++;

    for (int i = 0; i < thetas.size(); i++) {
      thetas[i] = abs(full_eqns[obj_row][i] / full_eqns[min_pos][i]);
      thetas_no_mod[i] = (full_eqns[obj_row][i] / full_eqns[min_pos][i]);
    }
    cout << "MIN_POSE = " << min_pos << endl;
    cout << "THETAS: ";
    print1D(thetas);
    thetas_cache.push_back(thetas_no_mod);
    int min_theta_pose = 0;

    for (int i = 0; i < thetas.size(); i++) {
      if (full_eqns[min_pos][i] > 0) continue;
      if (thetas[min_theta_pose] >= thetas[i]) {
        min_theta_pose = i;
      }
    }
    // cout<<"MIN ELE COL LOC"<<min_theta_pose<<endl;
    // cout<<"MIN ELE COL 2"<<thetas[min_theta_pose]<<endl;

    cout << "MIN_theta_POSE = " << min_theta_pose << endl;
    cout << "Pivot Element = " << full_eqns[min_pos][min_theta_pose] << endl;
    double pivot = full_eqns[min_pos][min_theta_pose];
    // break;
    // pair<int, int> pv(min_theta_pose, min_pos);
    pair<int, int> pv(min_pos, min_theta_pose);
    pivot_cache.push_back(pv);

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
  cout << "a: Number of Iters \nb: Leaving Vars \nc: Entering vars (deltaj/Xr) "
          "+ obj val\n"
          "d: Optimal Solution \nx: To end"
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
        cout << "Enter ith iter <= " << (count - 1) << endl;
        int i;
        cin >> i;
        // i = i + 1;
        // cout << "Non Basic Solutions" << endl;
        int position = pivot_cache[i].first;
        cout << "POSITION " << position << endl;
        print1D(index_cache[i].second);
        if (index_cache[i].second[position] <= index_cache[0].first.size())
          cout << "x" << index_cache[i].second[position] << "\t";
        else
          cout << "z" << (index_cache[i].second[position] -
                          index_cache[i].first.size())
               << "\t";
        cout << table_cache[i][position][full_eqns[0].size() - 1] << endl;
        break;
      }

      case 'c': {
        cout << "Enter ith iter <= " << (count - 1) << endl;
        int i;
        cin >> i;
        int position = pivot_cache[i].second;
        cout << "POSITION " << position << endl;
        print1D(index_cache[i].first);
        if (index_cache[i].first[position] <= index_cache[0].first.size())
          cout << "x" << index_cache[i].first[position] << "\t";
        else
          cout << "z"
               << (index_cache[i].first[position] - index_cache[i].first.size())
               << "\t";

        cout << "Min Ratio: " << thetas_cache[i][position] << endl;
        cout << "OBJECTIVE FUNCTION = " << table_cache[i][full_eqns.size() - 1]
                                                      [full_eqns[0].size() - 1]
             << endl;

        break;
      }

      case 'd': {
        if (unbounded) {
          cout << "Solution is Unbounded" << endl;
        }

        if (!unbounded) {
          cout << "Optimal Value: "
               << status *
                      full_eqns[full_eqns.size() - 1][full_eqns[0].size() - 1]
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
  if (!unbounded) {
    return status * full_eqns[full_eqns.size() - 1][full_eqns[0].size() - 1];
  } else {
    return -9000;
  }
}

int main(int argc, char const *argv[]) {
  cout << "Enter the number of Equations" << endl;
  int n, r;
  cin >> n;
  cout << "Enter the number of variables" << endl;
  cin >> r;
  vector<vector<double> > eqns;
  vector<int> var_record(r, 0);
  vector<int> signs;
  // vector<vector<double> > eqns_format(n, vector<double>((n + r + 1), 0.0));
  cout << "ENTER Variable types 1 for free 0 for not free" << endl;
  int count_num = 0;
  for (int i = 0; i < r; i++) {
    cin >> var_record[i];
    if (var_record[i] == 1) count_num++;
  }
  vector<vector<double> > eqns_only(n, vector<double>(r + 1 + count_num, 0.0));
  cout << "a1x1 + a2x2 + ... + arxr (>/</=)(specify when asked -1/+1/0) C"
       << endl;
  for (int i = 0; i < n; i++) {
    vector<double> arr(2 + r, 0);
    cout << "enter sign" << endl;
    cin >> arr[r];
    signs.push_back(arr[r]);
    int k = r;
    for (int j = 0; j < r; j++) {
      cin >> arr[j];
      // eqns_format[i][j] = arr[j];
      eqns_only[i][j] = arr[j] * arr[r];
      if (var_record[j] == 1) {
        eqns_only[i][k] = -arr[j] * arr[r];
        k++;
      }
    }

    // eqns_format[i][i + r] = arr[r];
    cout << "enter C" << endl;
    cin >> arr[r + 1];
    arr[r + 1] = arr[r + 1] * arr[r];
    // eqns_format[i][n + r] = arr[r + 1];
    eqns_only[i][r + count_num] = arr[r + 1];
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
    opt_full_splx[i] = -temp * eq_flag;
    optimal_eq.push_back(-temp * eq_flag);
  }
  for (int i = 0; i < r; i++) {
    if (var_record[i] == 1) optimal_eq.push_back(-opt_full_splx[i]);
  }

  vector<int> dual_var_record(n, 0);
  vector<int> dual_signs(r, 0);
  int dual_count = 0;

  for (int i = 0; i < n; i++) {
    if (signs[i] == 0) {
      dual_var_record[i] = 1;
      dual_count++;
    }
  }

  for (int i = 0; i < r; i++) {
    if (var_record[i] == 0) dual_signs[i] = 1;
  }
  vector<vector<double> > eqns_dual_only(
      r + count_num, vector<double>(n + 1 + dual_count, 0.0));
  int dk = r;
  for (int i = 0; i < r; i++) {
    int k = n;
    int k_dual = n;
    for (int j = 0; j < n; j++) {
      eqns_dual_only[i][j] = -eqns_only[j][i] * signs[j];
      if (dual_var_record[j] == 1) {
        eqns_dual_only[i][k] = -eqns_dual_only[i][j];
        k++;
      }
      if (var_record[i] == 1) {
        eqns_dual_only[dk][j] = eqns_only[j][i] * signs[j];
        if (dual_var_record[j] == 1) {
          eqns_dual_only[dk][k_dual] = -eqns_dual_only[i][j];
          k_dual++;
        }
      }
    }
    // cout<<"CURRENT I="<<i<<endl;
    eqns_dual_only[i][n + dual_count] = optimal_eq[i] * eq_flag;
    if (var_record[i] == 1) {
      eqns_dual_only[dk][n + dual_count] = -optimal_eq[i] * eq_flag;
      dk++;
    }
  }
  vector<double> optimal_dual_eq;
  for (int i = 0; i < n; i++) {
    optimal_dual_eq.push_back(eqns_only[i][eqns_only[0].size() - 1]);
  }
  optimal_dual_eq.push_back(0.0);
  optimal_eq.push_back(0.0);
  // cout<<"ENTER 1 for using primal equations, 2 for using dual
  // equations"<<endl;
  // int choice;
  // cin>>choice;
  double opt1 = 0.0, opt2 = 0.0;
  cout << "PRIMAL START" << endl;
  print2D(eqns_dual_only);
  print1D(opt_full_splx);
  print1D(optimal_eq);
  opt1 = simplex_solver(eqns_only, optimal_eq, eq_flag);
  cout << "PRIMAL END" << endl;
  cout << "DUAL START" << endl;
  print2D(eqns_dual_only);
  print1D(optimal_dual_eq);
  opt2 = simplex_solver(eqns_dual_only, optimal_dual_eq, -eq_flag);
  cout << "DUAL END" << endl;
  double max_out = (opt1 > opt2) ? opt1 : opt2;
  cout << "OPTIMAL VALUE OF THE EQUATION = " << max_out << endl;
  return 0;
}
