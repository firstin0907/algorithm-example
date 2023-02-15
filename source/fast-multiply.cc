#include <cstdio>
#include <cstring>
#include <cmath>

#include <vector>
#include <algorithm>
#include <complex>

using namespace std;

using cpx = complex <double>;
using vcpx = vector<cpx>;

constexpr int MAX_N = 1'111'111;
const double PI = acos(-1);

char a_buf[MAX_N], b_buf[MAX_N];

// transform char[] into vector<cpx> form.
void set_vector(const char* buf, int len, vcpx& v)
{
    int cnt = 0;
    for (int i = len - 1; i >= 0; i -= 2)
    {
        int val = buf[i] - '0';
        if (i != 0) val += (buf[i - 1] - '0') * 10;
        v[cnt++] = val;
    }
}

// conduct discrete fourier transform(DFT) or inverse DFT(IDFT).
// when rev = 1, f transforms as f(1), f(w), f(w^2), ..., f(w^n) where w = e^(2 * PI / n * i)
// when rev = -1, f transforms as f(1), f(w), f(w^2), ..., f(w^n) where w = e^(-2 * PI / n * i)
// where n is size of vector(degree of f)
void dft(vcpx& f, int rev)
{
    int n = f.size();

    // step 1. relocation the elements of f
    for (int i = 1; i < n; i++)
    {
        int opp = 0;
        for (int k = 1; k < n; k *= 2)
        {
            opp *= 2;
            if (i & k) ++opp;
        }
        if (i < opp) swap(f[i], f[opp]);
    }

    // step 2. get function value of f when applied every 1, w, w^2, ..., w^n
    for (int k = 2; k <= n; k *= 2)
    {
        double angle = rev * 2 * PI / k;
        cpx w(cos(angle), sin(angle));

        for (int s = 0; s < n; s += k)
        {
            cpx wp = 1;
            for (int i = s; i < s + k / 2; i++)
            {
                cpx even = f[i];
                cpx odd = f[i + k / 2];

                f[i] = even + odd * wp;
                f[i + k / 2] = even - odd * wp;
                wp *= w;
            }
        }
    }

    // step 3. if IDFT is being conducted, divide every elements by n, additionally.
    if(rev == -1)
    {
        for(int i = 0; i < n; i++) f[i] /= n;
    }
}

// multiply two vector<cpx> by using fast fourier transform(FFT).
// it consumes O(nlogn) time complexity.
void vcpx_multiply(vcpx& a, vcpx& b)
{
    const size_t n = a.size();

    dft(a, 1);
    dft(b, 1);

    for (int i = 0; i < n; i++) a[i] = a[i] * b[i];

    dft(a, -1);
}

void vcpx_to_array(vcpx& v, int arr[])
{
    const int n = v.size();

    for (int i = 0; i < n; i++) arr[i] = round(v[i].real());
    for (int i = 0; i < n - 1; i++) arr[i + 1] += arr[i] / 100, arr[i] %= 100;
}

// print result
void print_array(int arr[], int length)
{
    int flag = 0;
    for (int i = length - 1; i >= 0; i--)
    {
        if (flag) printf("%02d", arr[i]);
        else if (arr[i] != 0)
        {
            flag = 1;
            printf("%d", arr[i]);
        }
    }
    if (!flag) printf("0");
}

int answer[MAX_N << 1];
void fast_multiply(const char* a_buf, const char* b_buf)
{
    int alen = strlen(a_buf), blen = strlen(b_buf), n;
    for (n = 1; n <= alen / 2 + 1 || n <= blen / 2 + 1; n *= 2);
    n *= 2;

    vcpx a(n), b(n);
    set_vector(a_buf, alen, a);
    set_vector(b_buf, blen, b);

    vcpx_multiply(a, b);  
    
    vcpx_to_array(a, answer);
    print_array(answer, n);
}


int main()
{
    puts("input two non-negative large number (<= 10^100000) ");
    scanf("%s %s", a_buf, b_buf);
    fast_multiply(a_buf, b_buf);
    
    return 0;
}