#include <iostream>
#include <vector>

using namespace std;
vector<vector<int> > data;

void print2D(vector<vector<int> > mat)
{
    int n_rows = mat.size();
    int n_cols = mat[0].size();
    cout<<"N: "<<n_rows<<endl;
    for(int i=0; i<n_rows; i++)
    {
        for(int j=0; j<n_cols; j++)
            cout<<mat[i][j]<<"\t";
        cout<<endl;   
    }
}

void print1D(vector<int> vec)
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
        data.push_back(dat);
        return;
    }
    if (i >= n)
        return;
    dat[index] = arr[i];
    combinations(arr, dat, n, r, index, i+1);
    combinations(arr, dat, n, r, index+1, i+1);
}

int main(int argc, char const *argv[])
{
    vector<int> arr{1,2,3,4,5};
    int r = 3;
    int n = arr.size();
    vector<int> dat(r,0);
    combinations(arr, dat, n, r, 0, 0);
    cout<<"Computed combinations"<<endl;
    print2D(data);
    return 0;
}
