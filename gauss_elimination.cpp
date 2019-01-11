#include <iostream>
#include <cmath>
#include <vector>
#include <utility>

using namespace std;

void print2D(vector<vector<double> > mat)
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

void print1D(vector<double> vec)
{
    
    for(int i = 0; i < vec.size(); i++)
    {
        cout <<vec[i]<<" ";
    }
    cout<<endl;
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
            max.first = abs(max.first);
            for(int i = k+1; i < n_row; i++)
            {
                if(mat[i][l] != 0 && max.first<abs(mat[i][l]))
                {
                    max.first = abs(mat[i][l]);
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

int main(int argc, char const *argv[])
{
    vector<vector<double> > vect;
    cout<<"Enter the number of Equations"<<endl;
    int n, r;
    cin>>n;
    cout<<"Enter the number of variables"<<endl;
    cin>>r;
    cout<<"a1x1 + a2x2 + ... + arxr = C"<<endl;
    for(int i = 0; i < n; i++)
    {
        vector<double> arr(r+1, 0);
        for(int j = 0; j < r; j++)
            cin>>arr[j];
        cout<<"enter C"<<endl;
        cin>>arr[r]; 
        cout<<"Got Eqn: "<<i<<endl;
        vect.push_back(arr);
    }
    cout <<"Initial Mat"<<endl;
    print2D(vect);
    pair<vector<vector<double> >, vector<double> > p = gauss_elemination(vect);
    cout<<"Augmented Matrix"<<endl;
    print2D(p.first);
    cout<<"Solution Vector"<<endl;
    print1D(p.second);
    return 0;
}
