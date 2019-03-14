#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

using namespace std;
#define BigM 1e4
double global_best = 0;
int ele_num = 0;
vector<double> finale;
vector<vector<vector<double> > > final_tableu;
vector<pair<vector<int>, vector<int> > > final_indices;

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

pair<bool, pair<double, vector<double> > > simplex_solver(vector<vector<double> > full_eqns,
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
    if(full_eqns[i][cols - 2] == -1)
      M_index.push_back(-BigM);
    else
    {
      M_index.push_back(0.0);
    }
    
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
  //cout << "DONE" << endl;
  pair<bool, pair<double, vector<double> > > optim;
   bool infeasible = false;
   for (int j = 0; j < indices.second.size(); j++) {
          if (eqns_final[j][eqns_final[0].size() - 1] < 0) 
          infeasible = true;
        }

   if (infeasible || unbounded)
   {
       //cout<<"UBNOUNDED"<<endl;
       optim.first = false;
       return optim;
   }
   else
   {
       //cout<<"FEASIBLE"<<endl;
       optim.first = true;
   }

   if (optim.first)
   {
       optim.second.first = eqns_final[eqns_final.size() - 1][eqns_final[0].size() - 1];
       vector<double> v(ele_num, 0.0);
       cout<<"ENTERED"<<endl;
       for(int i = 0; i < indices.second.size(); i++)
       {
           //cout<<indices.second[i]<<","<<indices.first.size()<<endl;
           if (indices.second[i] <= indices.first.size())
           {
               v[indices.second[i]-1] = eqns_final[i][eqns_final[0].size() - 1];
           }
       }
       optim.second.second = v;
       cout<<"DONE"<<endl;
       //print2D(eqns_final);
       cout<<endl;
       final_tableu.push_back(eqns_final);
       final_indices.push_back(indices);
       //print1D(optim.second.second);
       //cout<<optim.second.first<<endl;
        return optim;
   }
  /*

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
  }*/
}

pair<vector<vector<double> >,pair<vector<int>, vector<int> > >  dual_simplex_solver(vector<vector<double> > full_eqns,
                      vector<double> full_obj, pair<vector<int>, vector<int> > p, int status) {
  // Initialization
  int n = full_eqns.size() - 1;         // rows
  int r = full_eqns[0].size() - 1;  // cols
  pair<vector<int>, vector<int> > indices = p;/*
  vector<int> sol_nb(r, 0);
  vector<int> sol_b(n, 0);
  for (int i = 0; i < sol_nb.size(); i++) sol_nb[i] = i + 1;  // the zero vars
  for (int i = 0; i < sol_b.size(); i++)
    sol_b[i] = i + 1 + r;  // the non zero vars
  indices = make_pair(sol_nb, sol_b);
  //full_eqns.push_back(full_obj);
  */
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
    cout<<"PP!"<<endl;
    for (int i = 0; i < full_eqns.size() - 1; i++) {
      if (full_eqns[min_pos][last_col] >= full_eqns[i][last_col] &&
          full_eqns[i][last_col] < 0) {
        min_pos = i;
        done = false;
      }
    } 
    cout<<"PP2"<<endl;
    unbounded = true;
    for (int i = 0; i < full_eqns[0].size(); i++) {
      if (full_eqns[min_pos][i] < 0) unbounded = false;
    }
    cout<<"PP3"<<endl;
    if (done || unbounded) {
      
        cout<<"VERY BAD"<<unbounded<<endl;
      if (done) unbounded = false;
      break;
    }

    else
      count++;
    cout<<"GOT YOU"<<endl;
    for (int i = 0; i < thetas.size(); i++) {
      thetas[i] = abs(full_eqns[obj_row][i] / full_eqns[min_pos][i]);
      thetas_no_mod[i] = (full_eqns[obj_row][i] / full_eqns[min_pos][i]);
    }
    cout << "MIN_POSE = " << min_pos << endl;
    cout << "THETAS: ";
    print1D(thetas);
    thetas_cache.push_back(thetas_no_mod);
    int min_theta_pose = 0;
    int min = 1e6;

    for (int i = 0; i < thetas.size(); i++) {
      if (full_eqns[min_pos][i] > 0) continue;
      if (min >= thetas[i]) {
        min_theta_pose = i;
        min = thetas[min_theta_pose];
      }
    }


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
        if (i == pv.first || j == pv.second) 
            continue;
        double p = 0.0, q = 0.0, r = 0.0, s = 0.0;
        p = pivot;
        s = full_eqns[i][j];
        r = full_eqns[pv.first][j];
        q = full_eqns[i][pv.second];
        full_eqns[i][j] = (p * s - q * r) / p;
        //cout<<full_eqns[i][j]<<" , "<<i<<" , "<<j<<endl;
      }
    }
    //print2D(full_eqns);
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
  //final_tableu = full_eqns;
  //final_indices = indices;
  cout<<"TABLE OF THE DAY"<<endl;
  //print_simplex_table(final_tableu, final_indices);
  return make_pair(full_eqns, indices);
  /*
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
  }*/
}

pair<bool, int> check_f_flag(vector<double> v, int k)
{
    vector<double> f(v.size(),0.0);
    double max = 0.0;
    int max_pose = 0;
    bool f_flag = true;
    for(int  i = k; i < v.size(); i++)
    {
      f[i] = v[i] - floor(v[i]);
      if (f[i]>0)
      {
        max_pose = i;
        break;
      }
    }
    
    for(int i = k; i < v.size(); i++)
    {
            f[i] = v[i] - floor(v[i]);
            if (f[i] != 0.0)
                f_flag = false;
            /*
            if (f[i]>max)
            {
                max = f[i];
                max_pose = i;
            }*/
    }
    return make_pair(f_flag, max_pose);
}
void branch_and_bound_recurse(vector<vector<double> > full_eqns,
                    vector<double> full_obj, int k)
{
    pair<bool, pair<double, vector<double> > > sol_init = simplex_solver(full_eqns, full_obj);
    vector<double> v = sol_init.second.second;
    if (k > full_eqns[0].size()-2)
        return;
    if (!sol_init.first)
    {
        cout<<"INFEASIBLE"<<endl;
        return;
    }
    pair<bool, int> status;
    status = check_f_flag(v, k);
    if (status.first)
    {
        cout<<"OPTIMAL of this iter"<<sol_init.second.first<<endl;
        if (sol_init.second.first > global_best)
        {
            global_best = sol_init.second.first;
            finale = sol_init.second.second;
        }
    }
    if(!status.first)
    {
        pair<double, double> limits;
        int pose = k;//status.second;
        limits.first = floor(v[pose]);
        limits.second = floor(v[pose]) + 1;
        vector<vector<double> > full_new1 = full_eqns;
        vector<vector<double> > full_new2 = full_eqns;
        vector<double> v1(full_eqns[0].size(), 0.0);
        v1[v1.size()-1] = limits.first;
        v1[pose] = 1;
        v1[v1.size()-2] = 1;
        full_new1.push_back(v1);
        vector<double> v2(full_eqns[0].size(), 0.0);
        v2[v2.size()-1] = limits.second;
        v2[pose] = 1;
        v2[v2.size()-2] = -1;
        full_new2.push_back(v2);
        k = pose + 1;
        branch_and_bound_recurse(full_new1, full_obj, k);
        branch_and_bound_recurse(full_new2, full_obj, k);
    }   
}

void cutting_plane(vector<vector<double> > full_eqns,
                    vector<double> full_obj)
{
  pair<bool, pair<double, vector<double> > > sol_init = simplex_solver(full_eqns, full_obj);
  cout<<"HAHAHAH"<<endl;
  print1D(sol_init.second.second);
  
  bool f_flag = true;
  int k = 0;
  //bool n_flag=true;
  vector<vector<double> > f_tableu;
  pair<vector<int>, vector<int> > f_indices;
  do
  {
    /* code */
  /*
  if (n_flag)
  {
    f_tableu = final_tableu;
    f_indices = final_indices;
    n_flag = false;
  }
  else
  {
    
  }*/
  
  vector<vector<double> > current_tableu;
  pair<vector<int>, vector<int> > current_indices;
  current_tableu = final_tableu[final_tableu.size() -1];
  current_indices = final_indices[final_indices.size() -1];
  vector<double> f(current_tableu.size() - 1,0.0);
  int col = current_tableu[0].size() - 1;
  int row = current_tableu.size() - 1;
  cout<<current_tableu[row][col]<<endl;
  k++;
  double max = 0.0;
  int max_pose = 0;
  f_flag = true;
  cout<<"F Diaries"<<endl;
  for(int i = 0; i < row; i++)
    {
      f[i] = current_tableu[i][col] - floor(current_tableu[i][col]);
      if (f[i] >= 1e-14 && f[i] < 0.9)
      {
        cout<<f[i]<<" , "<<i<<endl;
        f_flag = false;
      }
      if (f[i]>max)
        {
          max = f[i];
          max_pose = i;
        }
    }
    print1D(f);
    if(f_flag)
    {
      cout<<"DONE STATUS "<<f_flag<<endl;
      cout <<"FINALE"<<endl;
      print_simplex_table(current_tableu, current_indices);
      break;
    }
  vector<double> v(col + 1,0.0);
  for(int i = 0; i < v.size(); i++)
  {
    v[i] = -(current_tableu[max_pose][i] - floor(current_tableu[max_pose][i]));
  }
  current_tableu.insert(current_tableu.end()-1,v);
  cout<<"FOR SOME FUN"<<endl;
  print2D(current_tableu);
  //print2D(current_tableu);
  cout <<"THE FINAL"<<endl;
  current_indices.second.push_back(current_indices.first.size() + current_indices.second.size() + 1);
  cout<<"PRINTING TABLE"<<endl;
  print_simplex_table(current_tableu, current_indices);
  cout<<"STARTING DUAL"<<endl;
  pair<vector<vector<double> >,pair<vector<int>, vector<int> > > final_one = dual_simplex_solver(current_tableu, full_obj, current_indices, 1);
  final_tableu.push_back(final_one.first);
  final_indices.push_back(final_one.second);
  } while (!f_flag);
  
}
/*
void branch_and_bound(vector<vector<double> > full_eqns,
                    vector<double> full_obj)
{
    pair<bool, pair<double, vector<double> > > sol_init = simplex_solver(full_eqns, full_obj);
    
    vector<double> v = sol_init.second.second;
    
    if (!sol_init.first)
    {
        cout<<"INFEASIBLE"<<endl;
        return;
    }
        
    pair<bool, int> status;
    status = check_f_flag(v);
    while(!status.first)
    {
        pair<double, double> limits;
        int pose = status.second;
        limits.first = floor(v[pose]);
        limits.second = floor(v[pose]) + 1;
        vector<vector<double> > full_new1 = full_eqns;
        vector<vector<double> > full_new2 = full_eqns;
        vector<double> v1(full_eqns[0].size(), 0.0);
        v1[v1.size()-1] = limits.first;
        v1[pose] = 1;
        v1[v1.size()-2] = 1;
        full_new1.push_back(v1);
        vector<double> v2(full_eqns[0].size(), 0.0);
        v2[v2.size()-1] = limits.second;
        v2[pose] = 1;
        v2[v2.size()-2] = -1;
        full_new2.push_back(v2);
        cout<<"MAT1"<<endl;
        print2D(full_new1);
        pair<bool, pair<double, vector<double> > > sol_init1 = simplex_solver(full_new1, full_obj);
        cout<<"MAT2"<<endl;
        print2D(full_new2);
        pair<bool, pair<double, vector<double> > > sol_init2 = simplex_solver(full_new2, full_obj);

    }
    if (status.first)
    {
        cout<<"Optimal VAL: "<<sol_init.second.first<<endl;
        print1D(v);
        return;
    }
}*/

int main(int argc, char const *argv[]) {
  cout << "Enter the number of Equations" << endl;
  int n, r;
  cin >> n;
  cout << "Enter the number of variables" << endl;
  cin >> r;
  ele_num = r;
  cout << "a1x1 + a2x2 + ... + arxr (>/</=)(specify when asked -1/+1/0) C"
       << endl;
  vector<vector<double> > eqns;
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
  int num = 0;
  cout<<"Enter 1 for cutting plane \n2 for branch and bound method"<<endl;
  cin>>num;
  if (num == 2)
  {
    branch_and_bound_recurse(eqns, optimal_eq, 0);
     cout<<"THE BEST OF THE BEST"<<endl;
    cout<<"OPTIMAL VAL"<<global_best<<endl;
    print1D(finale);
  }
  if (num == 1)
    cutting_plane(eqns, optimal_eq);
  //simplex_solver(eqns, optimal_eq);
 
  return 0;
}
