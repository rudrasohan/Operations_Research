#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

using namespace std;

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

template <typename T>
vector<vector<T> > mat_prod(vector<vector<T> > A, vector<vector<T> > B)
{
    vector<vector<T> > M(A.size(), vector<T>(A.size(), 0.0));
    for(int i = 0; i < M.size(); i++)
    {
        for(int j = 0; j < M.size(); j++)
        {
            for(int k = 0; k < M.size(); k++)
            {
                M[i][j] += A[i][k] * B[k][j];
            }
            
        }
        
    }
    return M;
}

template <typename T>
vector<vector<T> > I(int n)
{
    vector<vector<T> > M(n, vector<T>(n, 0.0));
    for(int i = 0; i < n; i++)
        M[i][i] = 1;

    return M;
}


vector<double> eta(vector<double> e, int k)
{
    double pivot = e[k];
    for(int i = 0; i < e.size(); i++)
    {
        if(i == k)
            e[i] = 1/pivot;
        else
            e[i] /= -pivot;
        
    }
    return e;
}

template<typename T>
vector<T> vec_prod(vector<vector<T> > A, vector<T> v)
{
    vector<T> P(v.size(), 0.0);
    for(int i = 0; i < P.size(); i++)
    {
        for(int j = 0; j < P.size(); j++)
        {
            P[i] += A[i][j] * v[j];
        }
        
    }
    return P; 
}

template<typename T>
vector<T> extract_col(vector<vector<T> > A, int k)
{
    vector<T> col(A.size(), 0.0);
    for(int i = 0; i < A.size(); i++)
    {
        col[i] = A[i][k];   
    }
    return col;
}

template<typename T>
vector<vector<T> > subst_col(vector<vector<T> > A, vector<T> col, int k)
{
    for(int i = 0; i < col.size(); i++)
    {
        A[i][k] = col[i];   
    }
    return A;
}


vector<vector<double> > inv(vector<vector<double> > A)
{
    int n_iters = A.size();
    vector<vector<double> > B = I<double>(n_iters);
    vector<vector<double> > Id = I<double>(n_iters);
    vector<vector<double> > B_inv = B;
    for(int i = 0; i < n_iters; i++)
    {
        vector<double> C = extract_col(A, i);
        vector<double> e = vec_prod(B_inv, C);
        vector<double> et = eta(e, i);
        vector<vector<double> > E_i = subst_col(Id, et, i);
        B = subst_col(B, C, i);
        B_inv = mat_prod(E_i, B_inv);
    }
    return B_inv;
}

int main(int argc, char const *argv[])
{
    cout<<"Enter Mat Size()"<<endl;
    int n = 0;
    cin>>n;
    cout<<"Enter Elements of mat 1"<<endl;
    vector<vector<double> > A;
    vector<vector<double> > B;
    for(int i = 0; i < n; i++)
    {
        vector<double> temp(n, 0.0);
        for(int j = 0; j < n; j++)
        {
            cin>>temp[j];
        }
        A.push_back(temp);
    }
    print2D(inv(A));
    print2D(mat_prod(A, inv(A)));
    return 0;
}
