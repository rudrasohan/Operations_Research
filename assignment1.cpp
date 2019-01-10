#include <iostream>
#include <vector>
#include <utility>


using namespace std;

vector<vector<int> > comb_ind;

template <typename T>
void print2D(vector<vector<T> > mat)
{
    int n_rows = mat.size();
    int n_cols = mat[0].size();
    for(int i=0; i<n_rows; i++)
    {
        for(int j=0; j<n_cols; j++)
            cout<<mat[i][j]<<"\t";
        cout<<endl;   
    }
}

template <typename T>
void print1D(vector<T> vec)
{
    
    for(int i = 0; i < vec.size(); i++)
    {
        cout <<vec[i]<<" ";
    }
    cout<<endl;
}

void combinations(vector<int> arr, vector<int> dat, int n, int r, int index, int i)
{
    if (index == r)
    {
        comb_ind.push_back(dat);
        return;
    }
    if (i >= n)
        return;
    dat[index] = arr[i];
    combinations(arr, dat, n, r, index, i+1);
    combinations(arr, dat, n, r, index+1, i+1);
}

vector<vector<double> > sub(vector<vector<double> > mat, double factor, int next, int k)
{
    int n_cols = mat[0].size();
    
    for(int  i = 0; i < n_cols; i++)
    {
        mat[next][i] += factor * mat[k][i];
    }
    return mat;
}

vector<vector<double> > swap_row(vector<vector<double> > mat, int i, int j)
{
    int N = mat[0].size();
    for(int k=0; k<N; k++)
    {
        double temp = 0.0;
        temp = mat[i][k];
        mat[i][k] = mat[j][k];
        mat[j][k] = temp;
    }
    return mat;
}

vector<double> back_substitution(vector<vector<double> > mat)
{
    int n_row = mat.size();
    int n_col = mat[0].size();
    cout<<n_row<<","<<n_col<<endl;
    vector<double> x(n_row, 0.0);
    for(int i = n_row-1; i >=0 ; i--)
    {
        double cache = 0.0;
        for(int j = i+1; j < n_col-1; j++)
        {
            cache += mat[i][j]*x[j];
        }
        x[i] = (mat[i][n_col-1] - cache)/mat[i][i];
    }
    return x;
    
}

pair<vector<vector<double> >, vector<double> > gauss_elemination(vector<vector<double> > mat)
{
    int n_row = mat.size();
    int n_col = mat[0].size();
    int k = 0, l = 0;
    while(l < n_col-2)
    {
        if(mat[k][l] != 0)
        {
            
            for(int i = k+1; i < n_row; i++)
            {
                if (mat[i][l] != 0)
                {
                    double factor = - mat[i][l] / mat[k][l];
                    mat = sub(mat, factor, i, k);
                    //print2D(mat);
                }
            }
            k++;
            l++;
        }
        else
        {
            bool column_change = true;
            pair<double,int> max (mat[k][l], k);
            for(int i = k+1; i < n_row; i++)
            {
                if(mat[i][l] != 0 && max.first<mat[i][l])
                {
                    max.first = mat[i][l];
                    max.second = i;
                    column_change = false;
                }
            }
            if (k != max.second && !column_change)
                mat = swap_row(mat, k, max.second);

            if (column_change)
                l++;
        }    
    }
    return make_pair(mat, back_substitution(mat));
}

vector<vector<double> > create_data_aug_mat(vector<vector<double> > full_mat, vector<int> indices)
{
    int M = full_mat.size();
    int N = full_mat[0].size();
    vector<vector<double> > mat(M, vector<double>(M+1, 0.0));
    for(int i=0; i<indices.size(); i++)
    {
        if(indices[i] < (N-2))
        {
            for(int j = 0; j < M; j++)
                mat[j][i] = full_mat[j][indices[i]];
        }
        else
        {
            int col_pos = N-2;
            int row_pos = indices[i] - (N-2);
            mat[i][i] = full_mat[row_pos][col_pos];
        }
    }
    
    for(int i = 0; i < M; i++)
    {
        mat[i][M] = full_mat[i][N-1];
    }
    return mat;
}

int main(int argc, char const *argv[])
{
    cout<<"Enter the number of Equations"<<endl;
    int n, r;
    cin>>n;
    cout<<"Enter the number of variables"<<endl;
    cin>>r;
    cout<<"a1x1 + a2x2 + ... + arxr (>/</=)(specify when asked -1/+1/0) C"<<endl;
    vector<vector<double> > eqns;
    for(int i = 0; i < n; i++)
    {
        vector<double> arr(r+2, 0);
        for(int j = 0; j < r; j++)
            cin>>arr[j];
        cout<<"enter sign"<<endl;
        cin>>arr[r];
        cout<<"enter C"<<endl;
        cin>>arr[r+1]; 
        eqns.push_back(arr);
    }
    cout<<"Equations Data"<<endl;
    print2D<double>(eqns);

    cout<<"Creating Combinations table"<<endl;
    int n_var = r + n;
    int m_val = n;
    vector<int> v(n_var, 0);
    for(int i = 0; i < n_var; i++)
        v[i] = i;
    vector<int> dat(m_val, 0);
    combinations(v, dat, n_var, m_val, 0, 0);
    print2D<int>(comb_ind);
    int comb = comb_ind.size();

    vector<vector<double> > basic_sols;
    for(int i=0; i<comb; i++)
    {
        vector<int> indices = comb_ind[i];
        vector<vector<double> > mat_test = create_data_aug_mat(eqns, indices);
        cout<<"PRINT TEST MAT"<<endl;
        print2D<double>(mat_test);
    }
    return 0;
}
