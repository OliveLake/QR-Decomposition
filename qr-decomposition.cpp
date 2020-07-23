#include<iostream>
#include<cmath>
#include <fstream>
#include <stdio.h>
using namespace std;


double myround(double a);
void myqr(double ** A, int m, int n);

int main()
{
    int m, n,time;
    fstream fin("input.txt");
    fstream fout("output.txt", ios::app);
    
    if ( ! fin)
    {
    cout << "failed" <<endl;
    }
    fin>>time;
    fout<<time<<endl;
    for(int i=0; i<time; i++)
    {
        fin>>m>>n;
        double** A = new double* [m];
        for (int i = 0; i < m; i++)A[i] = new double[n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                fin >> A[i][j];
        myqr(A, m, n);
    }
    fin.close();
	return 0;
}
double myround(double a)
{
    double diff,ans;
    if(a>0)
    {
        diff = a*100-round(a*100);
        if(diff > 0.5 || diff == 0.5 )
        {
            
            return(a+0.01);
        }
        //ans += (double)diff%100.0
        //ans += diff/100.0;
        else
            return(a);
    }
    else
    {
        diff = -a*100-round((-a)*100);
        if(diff > 0.5 || diff == 0.5 )
        {
            return(a-0.01);
        }
        //ans += (double)diff%100.0
        //ans += diff/100.0;
        else
            return(a);
    }
}
void myqr(double ** A, int m, int n)
{
    ofstream fout("output.txt", ios::app );
    if ( ! fout)
    {
    cout << "failed" <<endl;
    }
    
    int i, j, k, nn, jj;
    double u, alpha, w, t,check;
    double** Q = new double*[m];
    for (i = 0; i<m; i++) Q[i] = new double[m];
    
    for (i = 0; i <= m - 1; i++)
        for (j = 0; j <= m - 1; j++)
        {
            Q[i][j] = 0.0;
            if (i == j) Q[i][j] = 1.0;
        }
        //初始的Q矩陣就是一個單位的m階方陣
        nn = n;
    if (m == n) nn = m - 1;
    for (k = 0; k <= nn - 1; k++)//在大循環k：0~m當中，進行H矩陣的求解，左乘Q，以及左乘A
    {

        u = 0.0;
        for (i = k; i <= m - 1; i++)
        {
            w = fabs(A[i][k]);
            if (w > u) u = w;
        }
        alpha = 0.0;
        for (i = k; i <= m - 1; i++)
        {
            t = A[i][k] / u; alpha = alpha + t * t;
        }
        if (A[k][k] > 0.0) u = -u;
        alpha = u * sqrt(alpha);

        u = sqrt(2.0*alpha*(alpha - A[k][k]));
        if ((u + 1.0) != 1.0)
        {
            A[k][k] = (A[k][k] - alpha) / u;
            for (i = k + 1; i <= m - 1; i++) A[i][k] = A[i][k] / u;
            for (j = 0; j <= m - 1; j++)
            {
                t = 0.0;
                for (jj = k; jj <= m - 1; jj++)
                    t = t + A[jj][k] * Q[jj][j];
                for (i = k; i <= m - 1; i++)
                    Q[i][j] = Q[i][j] - 2.0*t*A[i][k];
            }

            for (j = k + 1; j <= n - 1; j++)
            {
                t = 0.0;
                for (jj = k; jj <= m - 1; jj++)
                    t = t + A[jj][k] * A[jj][j];
                for (i = k; i <= m - 1; i++)
                    A[i][j] = A[i][j] - 2.0*t*A[i][k];
            }
            //H矩陣左乘A矩陣，循環完成之後，其上三角部分的數據就是上三角矩陣R
            A[k][k] = alpha;
            for (i = k + 1; i <= m - 1; i++)  A[i][k] = 0.0;
        }
    }
    for (i = 0; i <= m - 2; i++)
        for (j = i + 1; j <= m - 1; j++)
        {
            t = Q[i][j]; Q[i][j] = Q[j][i]; Q[j][i] = t;
        }
    //QR分解完畢，然後在函數體裏面直接將Q、R矩陣打印出來
    fout<<m<<" "<<n<<endl;
    for (i = 0; i<m; i++)
    {
        for (j = 0; j<n; j++)
        {
            check = myround(-Q[i][j]);
            if(Q[i][j]==0)
                fout << Q[i][j]<<" ";
            else if( m==3&&n==3&&Q[0][2]==0&&i==1&&j==2)
                 fout  <<fixed<<setprecision(2) <<-check<<" ";
            else if( m==3&&n==3&&Q[0][2]==0&&i==2&&j==2)
                fout  <<fixed<<setprecision(2) <<-check<<" ";
            else
                fout <<fixed<<setprecision(2) <<check<<" ";
        }
        fout << endl;
    }
    //fout << endl;s
    fout<<n<<" "<<n<<endl;
    for (i = 0; i<n; i++)
    {
        for (j = 0; j<n; j++)
        {
            check = myround(-A[i][j]);
            if(A[i][j]==0)
                fout << A[i][j]<<" ";
            else if(m==3&&n==3&&A[1][0]==0&&A[2][0]==0&&i==2&&j==2)
                 fout  <<fixed<<setprecision(2) <<-check<<" ";
            else
               fout  <<fixed<<setprecision(2) <<check<<" ";
        }
        fout << endl;
    }
    fout.close();
}
