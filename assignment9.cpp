#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

using namespace std;

vector<vector<int> > comb_ind;
bool stable;
double most_min;

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

vector<vector<double> > simplex_solver(vector<vector<double> > full_eqns,
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
          if (thetas_cache[i - 1][i] < 0) continue;
          if (thetas_cache[i - 1][min_theta_pose] >= thetas_cache[i - 1][i]) {
            min_theta_pose = i;
          }
        }
        cout << thetas_cache[i - 1][min_theta_pose] << endl;
        break;
      }

      case 'd': {
        cout << "Enter i for ith table" << endl;
        int i;
        cin >> i;
        print_simplex_table(table_cache[i], index_cache[i]);
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
          double v =
              1.0 / full_eqns[full_eqns.size() - 1][full_eqns[0].size() - 1];
          if (most_min > 0)
            cout << "Optimal Value: " << v << endl;
          else
            cout << "Optimal Value: " << v + most_min << endl;

          vector<double> sols(indices.first.size(), 0.0);
          for (int i = 0; i < indices.second.size(); i++) {
            if (indices.second[i] <= indices.first.size()) {
              sols[indices.second[i] - 1] =
                  full_eqns[i][full_eqns[0].size() - 1] * v;
            }
          }
          cout << "OPTIMAL B POLICY" << endl;
          print1D(sols);
        }
        break;
      }

      default: { cout << "Invalid Option"; }
    }
  }
  return full_eqns;
}

double dual_simplex_solver(vector<vector<double> > full_eqns,
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
          double v = 1.0 /
                     full_eqns[full_eqns.size() - 1][full_eqns[0].size() - 1] *
                     status;
          if (most_min > 0)
            cout << "Optimal Value: " << v << endl;
          else
            cout << "Optimal Value: " << v + most_min << endl;

          vector<double> sols(indices.first.size(), 0.0);
          for (int i = 0; i < indices.second.size(); i++) {
            if (indices.second[i] <= indices.first.size()) {
              sols[indices.second[i] - 1] =
                  full_eqns[i][full_eqns[0].size() - 1] * v;
            }
          }
          cout << "OPTIMAL A POLICY" << endl;
          print1D(sols);
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

double get_max_min(vector<double> arr, int val) {
  int len = arr.size();
  int var = arr[0];
  for (int i = 0; i < len; i++) {
    if (val * arr[i] <= val * var) {
      var = arr[i];
    }
  }
  return var;
}

vector<vector<double> > add_2table(vector<vector<double> > mat, double k) {
  for (int i = 0; i < mat.size(); i++) {
    for (int j = 0; j < mat[0].size(); j++) {
      mat[i][j] += k;
    }
  }
  return mat;
}

vector<vector<double> > check_table(vector<vector<double> > mat) {
  int row, col;
  row = mat.size();
  col = mat[0].size();
  vector<double> B_loss(row, 0.0);
  vector<double> A_gain(col, 0.0);
  double seta = 0.0, setb = 0.0;
  for (int i = 0; i < row; i++) {
    B_loss[i] = get_max_min(mat[i], 1);
  }
  setb = get_max_min(B_loss, -1);
  print1D(B_loss);
  most_min = get_max_min(B_loss, 1);
  cout << "MOST_MIN: " << most_min << endl;
  for (int i = 0; i < col; i++) {
    vector<double> temp_array;
    for (int j = 0; j < row; j++) {
      temp_array.push_back(mat[j][i]);
    }

    A_gain[i] = get_max_min(temp_array, -1);
  }
  seta = get_max_min(A_gain, 1);
  print1D(A_gain);
  if (seta == setb) {
    stable = true;
    cout << "STABLE" << endl;
    cout << "OPTIMAL VALUE: " << stable << endl;
  } else {
    stable = false;
    cout << "UNSTABLE" << endl;
  }
  vector<vector<double> > mat_new;
  if (most_min < 0) {
    mat = add_2table(mat, -most_min);
  }
  return mat;
}

vector<vector<double> > create_eqns(vector<vector<double> > mat, int dual) {
  vector<vector<double> > tables(
      mat.size(), vector<double>(mat[0].size() + 1, (1.0 * dual)));
  for (int i = 0; i < mat.size(); i++) {
    for (int j = 0; j < mat[0].size(); j++) {
      tables[i][j] = mat[i][j] * dual;
    }
  }
  return tables;
}

vector<vector<double> > get_transpose(vector<vector<double> > mat) {
  vector<vector<double> > mat_T(mat[0].size(), vector<double>(mat.size(), 0.0));
  for (int i = 0; i < mat.size(); i++) {
    for (int j = 0; j < mat[0].size(); j++) {
      mat_T[j][i] = mat[i][j];
    }
  }
  return mat_T;
}

int main(int argc, char const *argv[]) {
  int n, r;
  cout << "Enter no of strategies for A" << endl;
  cin >> n;
  cout << "Enter no of strategies for B" << endl;
  cin >> r;
  cout << "ENTER DATA row-wise" << endl;
  vector<vector<double> > eqns_only(n, vector<double>(r, 0.0));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < r; j++) {
      cin >> eqns_only[i][j];
    }
    cout << "Got Row " << (i + 1) << endl;
  }

  cout << "Table Data" << endl;
  print2D<double>(eqns_only);
  vector<vector<double> > eqn = check_table(eqns_only);
  print2D(eqn);
  if (!stable) {
    vector<vector<double> > table = create_eqns(eqn, 1);
    print2D(table);
    vector<double> optimal_eq;
    // vector<double> opt_full_splx(n + r, 0.0);
    for (int i = 0; i < r; i++) {
      double temp = 0.0;
      // cin >> temp;
      // opt_full_splx[i] = -temp;
      optimal_eq.push_back(-1.0);
    }
    optimal_eq.push_back(0.0);
    // print1D(opt_full_splx);
    print1D(optimal_eq);
    vector<vector<double> > Btable = simplex_solver(table, optimal_eq);
    vector<vector<double> > Ttable = get_transpose(eqn);
    print2D(Ttable);
    vector<vector<double> > dual_eqn = create_eqns(Ttable, -1);
    print2D(dual_eqn);
    vector<double> dual_optimal_eq;
    for (int i = 0; i < n; i++) {
      double temp = 0.0;
      // cin >> temp;
      // opt_full_splx[i] = -temp;
      dual_optimal_eq.push_back(1.0);
    }
    dual_optimal_eq.push_back(0.0);
    dual_simplex_solver(dual_eqn, dual_optimal_eq, -1);
  }
  return 0;
}
