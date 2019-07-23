# include <iostream>
using namespace std;

//-----------------------コピペここから--------------------------

/* 拡張ユークリッドの互除法(ax + by = 1 を満たすx, y) */
template<typename T>
T extgcd(T a, T b, T& x, T& y) {
	T d = a;
	if (b != 0) {
		d = extgcd(b, a % b, y, x);
		y -= (a / b) * x;
	}
	else {
		x = 1; y = 0;
	}
	return d;
}

/* MODを法とする環 */
template<typename T, T MOD = 1000000007>
struct Mint {

	static constexpr T mod = MOD;

	T v;

	/* コンストラクタ群 */
	Mint() :v(0) {}

	Mint(signed _v) :v(_v) {}

	Mint(long long t) {
		v = t % mod;
		if (v < 0) v += mod;
	}

	/* べき乗を計算 */
	Mint pow(long long k) {
		Mint res(1), tmp(v);
		while (k) {
			if (k & 1) res *= tmp;
			tmp *= tmp;
			k >>= 1;
		}
		return res;
	}

	/* 単位元 */
	static Mint add_identity() { return Mint(0); }
	static Mint mul_identyty() { return Mint(1); }

	/* 逆元を求める */
	// オイラーの定理
	Mint inv() { return pow(mod - 2); }

	// 拡張ユークリッド互除法
	Mint inv_euclid() {
		T x, y;
		extgcd(v, mod, x, y);
		return Mint(x);
	}

	/* 演算子群 */
	Mint& operator+=(Mint a) {
		v += a.v;
		if (v >= mod) v -= mod;
		return *this;
	}

	Mint& operator-=(Mint a) {
		v += mod - a.v;
		if (v >= mod) v -= mod;
		return *this;
	}

	Mint& operator*=(Mint a) {
		v = 1LL * v * a.v % mod;
		return *this;
	}

	Mint& operator/=(Mint a) {
		return (*this) *= a.inv();
		// ユークリッドの互除法を使う場合
		// return (*this) *= a.inv_euclid();
	}

	Mint& operator+=(T a) { return (*this) += Mint(a); }
	Mint& operator-=(T a) { return (*this) -= Mint(a); }
	Mint& operator*=(T a) { return (*this) *= Mint(a); }
	Mint& operator/=(T a) { return (*this) /= Mint(a); }

	Mint operator+(Mint a) const { return Mint(v) += a; }
	Mint operator-(Mint a) const { return Mint(v) -= a; }
	Mint operator*(Mint a) const { return Mint(v) *= a; }
	Mint operator/(Mint a) const { return Mint(v) /= a; }

	Mint operator+(T a) const { return Mint(v) += a; }
	Mint operator-(T a) const { return Mint(v) -= a; }
	Mint operator*(T a) const { return Mint(v) *= a; }
	Mint operator/(T a) const { return Mint(v) /= a; }

	Mint operator-() const { return v ? Mint(mod - v) : Mint(v); }

	bool operator==(const Mint a) const { return v == a.v; }

	bool operator!=(const Mint a) const { return v != a.v; }

	bool operator<(const Mint a) const { return v < a.v; }

	bool operator>(const Mint a) const { return v > a.v; }

	/* 組み合わせ */
	static Mint comb(long long n, int k) {
		Mint num(1), dom(1);
		for (int i = 0; i < k; i++) {
			num *= Mint(n - i);
			dom *= Mint(i + 1);
		}
		return num / dom;
	}

};

template<typename T, T MOD> constexpr T Mint<T, MOD>::mod;

/* ostreamで出力できるように */
template<typename T, T MOD>
ostream& operator<<(ostream &os, Mint<T, MOD> m) { os << m.v; return os; }

//-----------------------コピペここまで--------------------------
