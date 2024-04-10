#include <iostream>

using namespace std;
void Comb(int* a, int n, int k);


void GenerateCases()
{
    int n = 12,k = 8;
    int *a = new int[n]();
    for (int i = 0; i < n; i++)
    {
     a[i] = i+1;
    }
    Comb(a, n, k);
    delete[] a;
}

void printRes(int* a, bool* index, int n)
{
    for (int i=0;i<n;i++)
    {
        if (index[i])
        {
            cout << a[i] << " ";
        }
    }
    cout << endl;
}

//检查最后k个位置是否已全变成0
bool hasDone(bool* index, int n, int k)
{
    for (int i=n-1;i>=n-k;i--)
    {
        if (!index[i])
        {
            return false;
        }
    }
    return true;
}

void Comb(int* a, int n, int k)
{
    bool *index = new bool[n]();
    //选中前k个位置
    for (int i = 0; i < k; i++)
    {
        index[i] = true;
    }
    printRes(a, index, n);
while (!hasDone(index, n, k))
    {
        //从左到右扫描数组
        for (int i = 0; i < n - 1; i++)
        {
            //找到第一个“10”组合将其变成"01"组合
            if (index[i] && !index[i + 1])
            {
                index[i] = false;
                index[i + 1] = true;

                //将"01"组合左边的1移到最左边
                int count = 0;
                for (int j = 0; j < i; j++)
                {
                    if (index[j])
                    {
                        index[j] = false;
                        index[count++] = true;
                    }
                }
                printRes(a, index, n);
                break;
            }
        }
    }
    delete[] index;
}


