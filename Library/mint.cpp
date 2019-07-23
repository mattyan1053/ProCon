using namespace std;

template<typename T, T MOD = 1000000007>
struct Mint {

    static constexpr T mod = MOD;

    T v;
    
    /* コンストラクタ群 */
    Mint():v(0){}

    Mint(signed _v):v(_v){}

    Mint(long long t){
        v = t % mod;
        if (v < 0) v += mod;
    }

    Mint pow(long long k){
        Mint res(1), tmp(v);
        while(k) {
            if(k & 1) res *= tmp;
            tmp *= tmp;
            k >> 1;
        }
        return res;
    }

    static Mint add_identity(){return Mint(0);}
    static Mint mul_identyty(){return Mint(1);}
    
    /* 逆元を求める */
    // 拡張ユークリッド互除法
    Mint inv(){
        Mint u(1), v(0), tmp(v), b(mod);
        while(b.v){
            Mint t = tmp / b;
            tmp -= t * b;
            u -= t * v;
            swap(u, v);
        }
        return u;
    }
    // オイラーの定理
    Mint inv_euler(){return pow(mod - 2);}

    /* 演算子群 */
    Mint& operator+=(Mint a){
        v += a.v;
        if (v >= mod) v -= mod;
        return *this;
    }

    Mint& operator-=(Mint a){
        v += mod - a.v;
        if (v >= mod) v -= mod;
        return *this;
    }

    Mint& operator*=(Mint a){
        v = 1LL * v * a.v % mod;
        return *this;
    }

    Mint& operator/=(Mint a){
        return (*this) *= a.inv();
    }

    Mint operator+(Mint a) const {return Mint(v) += a;}

    Mint operator-(Mint a) const {return Mint(v) -= a;}

    Mint operator*(Mint a) const {return Mint(v) *= a;}

    Mint operator/(Mint a) const {return Mint(v) /= a;}

    Mint operator-() const {return v ? Mint(mod - v) : Mint(v);}

    bool operator==(const Mint a) const {return v == a.v;}

    bool operator!=(const Mint a) const {return v != a.v;}

    bool operator<(const Mint a) const {return v < a.v;}

    bool operator>(const Mint a) const {return v > a.v;}

    /* 組み合わせ */
    static Mint comb(long long n, int k) {
        Mint num(1), dom(1);
        for(int i=0;i < k; i++){
            num *= Mint(n - i);
            dom *= Mint(i + 1);
        }
        return num / dom;
    }

};

template<typename T, T MOD> constexpr T Mint<T, MOD>::mod;
template<typename T, T MOD>
ostream& operator<<(ostream &os, Mint<T, MOD> m){ os << m.v; return os;}
