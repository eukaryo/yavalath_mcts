
#include<iostream>
#include<vector>
#include<string>
#include<map>
#include<set>
#include<queue>
#include<algorithm>
#include<cassert>
#include<cstdint>
#include<regex>
#include<random>
#include<cstdint>
#include<iomanip>
#include<thread>
#include<chrono>
#include<array>
#include<bitset>
#include<functional>
#include<exception>

//#define ONLY_BASIC_INSTRUCTION

#ifndef ONLY_BASIC_INSTRUCTION
#include <mmintrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include <tmmintrin.h>
#include <smmintrin.h>
#include <nmmintrin.h>
#include <wmmintrin.h>
#include <immintrin.h>
#endif

inline uint64_t popcount64_plain(uint64_t x) noexcept {
	x = x - ((x >> 1) & 0x5555'5555'5555'5555ULL);
	x = (x & 0x3333'3333'3333'3333ULL) + ((x >> 2) & 0x3333'3333'3333'3333ULL);
	x = (x + (x >> 4)) & 0x0F0F'0F0F'0F0F'0F0FULL;
	return (x * 0x0101'0101'0101'0101ULL) >> 56;
}

inline uint64_t popcount64(uint64_t x) noexcept {

#ifdef ONLY_BASIC_INSTRUCTION
	return popcount64_plain(x);
#else
	return _mm_popcnt_u64(x);
#endif

}

inline uint32_t bitscan_forward64(unsigned long *dest, const uint64_t x) noexcept {

	//xが非ゼロなら、立っているビットのうち最下位のものの位置をdestに代入して、非ゼロの値を返す。
	//xがゼロなら、ゼロを返す。このときのdestの値は未定義である。

#ifdef ONLY_BASIC_INSTRUCTION
	if (x) {
		*dest = popcount64(~x & (x - 1));
		return 1;
	}
	else {
		return 0;
	}
#elif _MSC_VER
	return _BitScanForward64(dest, x);
#else
	if (x) {
		*dest = __builtin_ctzll(x);
		return 1;
	}
	else {
		return 0;
	}
#endif

}


/*
void prepare_pdep_pext_tables() {
	uint8_t _pdep_table[16][16] = {};
	uint8_t _pext_table[16][16] = {};
	uint8_t _popcount_table[16] = {};

	const auto _pdep_naive = [](uint64_t a, uint64_t mask) {

		uint64_t answer = 0, k = 0;

		for (uint64_t m = 0; m < 64; ++m) {
			if (mask & (1ULL << m)) {
				if (a & (1ULL << k++)) {
					answer += (1ULL << m);
				}
			}
		}

		return answer;
	};

	const auto _pext_naive = [](uint64_t a, uint64_t mask) {

		uint64_t answer = 0, k = 0;

		for (uint64_t m = 0; m < 64; ++m) {
			if (mask & (1ULL << m)) {
				if (a & (1ULL << m)) {
					answer += (1ULL << k);
				}
				++k;
			}
		}

		return answer;
	};

	const auto popcount_naive = [](uint64_t x) {

		uint64_t answer = 0;

		for (uint64_t i = 0; i < 64; ++i) {
			if (x & (1ULL << i))++answer;
		}

		return answer;
	};

	const auto _init_tables = [&]() {
		for (uint64_t x = 0; x < 16; ++x)for (uint64_t y = 0; y < 16; ++y) {
			const uint64_t p = _pdep_naive(x, y);
			assert(p < 16ULL);
			_pdep_table[x][y] = (uint8_t)p;
		}

		for (uint64_t x = 0; x < 16; ++x)for (uint64_t y = 0; y < 16; ++y) {
			const uint64_t p = _pext_naive(x, y);
			assert(p < 16ULL);
			_pext_table[x][y] = (uint8_t)p;
		}

		for (uint64_t x = 0; x < 16; ++x) {
			const uint64_t p = popcount_naive(x);
			assert(p <= 4ULL);
			_popcount_table[x] = (uint8_t)p;
		}
	};

	const auto _pext_plain = [&](uint64_t a, uint64_t mask) {

		uint64_t answer = 0;

		for (uint64_t m = 0; mask; mask >>= 4, a >>= 4) {
			const uint64_t b = mask % 16;
			answer += ((uint64_t)_pext_table[a % 16][b]) << m;
			m += _popcount_table[b];
		}

		return answer;
	};

	const auto _pdep_plain = [&](uint64_t a, uint64_t mask) {

		uint64_t answer = 0;

		for (uint64_t m = 0; mask; mask >>= 4, m += 4) {
			const uint64_t b = mask % 16;
			answer += ((uint64_t)_pdep_table[a % 16][b]) << m;

			a >>= _popcount_table[b];
		}

		return answer;
	};

	const auto _test_pdep_pext = [&](uint64_t seed, int rep) {

		mt19937_64 rnd(seed);

		for (int i = 0; i < rep; ++i) {
			const uint64_t a = rnd(), mask = rnd();

			const uint64_t d0 = _pdep_naive(a, mask);
			const uint64_t d1 = _pdep_plain(a, mask);

			const uint64_t e0 = _pext_naive(a, mask);
			const uint64_t e1 = _pext_plain(a, mask);

			if (d0 != d1) {
				std::cout << "test failed!" << std::endl;
				std::cout << std::bitset<64>(d0) << std::endl;
				std::cout << std::bitset<64>(d1) << std::endl;
				break;
			}
			if (e0 != e1) {
				std::cout << "test failed!" << std::endl;
				std::cout << std::bitset<64>(e0) << std::endl;
				std::cout << std::bitset<64>(e1) << std::endl;
				break;
			}
		}
	};

	const auto _output_tables = [&]() {

		std::cout << "alignas(64) uint8_t pdep_table[16][16] = {" << std::endl;
		for (int i = 0; i < 16; ++i) {
			std::cout << "{";
			for (int j = 0; j < 16; ++j) {
				std::cout << int(_pdep_table[i][j]) << (j != 15 ? ", " : "}");
			}
			if (i != 15)std::cout << ",";
			std::cout << std::endl;
		}
		std::cout << "};" << std::endl;

		std::cout << "alignas(64) uint8_t pext_table[16][16] = {" << std::endl;
		for (int i = 0; i < 16; ++i) {
			std::cout << "{";
			for (int j = 0; j < 16; ++j) {
				std::cout << int(_pext_table[i][j]) << (j != 15 ? ", " : "}");
			}
			if (i != 15)std::cout << ",";
			std::cout << std::endl;
		}
		std::cout << "};" << std::endl;

		std::cout << "alignas(64) uint8_t popcount_table[16] = {";
		for (int i = 0; i < 16; ++i) {
			std::cout << int(_popcount_table[i]) << (i != 15 ? ", " : "};");
		}
		std::cout << std::endl;
	};

	_init_tables();
	_test_pdep_pext(12345, 10000);
	_output_tables();
}
*/


alignas(64) uint8_t pdep_table[16][16] = {
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 1, 2, 1, 4, 1, 2, 1, 8, 1, 2, 1, 4, 1, 2, 1},
{0, 0, 0, 2, 0, 4, 4, 2, 0, 8, 8, 2, 8, 4, 4, 2},
{0, 1, 2, 3, 4, 5, 6, 3, 8, 9, 10, 3, 12, 5, 6, 3},
{0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 8, 0, 8, 8, 4},
{0, 1, 2, 1, 4, 1, 2, 5, 8, 1, 2, 9, 4, 9, 10, 5},
{0, 0, 0, 2, 0, 4, 4, 6, 0, 8, 8, 10, 8, 12, 12, 6},
{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 7},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8},
{0, 1, 2, 1, 4, 1, 2, 1, 8, 1, 2, 1, 4, 1, 2, 9},
{0, 0, 0, 2, 0, 4, 4, 2, 0, 8, 8, 2, 8, 4, 4, 10},
{0, 1, 2, 3, 4, 5, 6, 3, 8, 9, 10, 3, 12, 5, 6, 11},
{0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 8, 0, 8, 8, 12},
{0, 1, 2, 1, 4, 1, 2, 5, 8, 1, 2, 9, 4, 9, 10, 13},
{0, 0, 0, 2, 0, 4, 4, 6, 0, 8, 8, 10, 8, 12, 12, 14},
{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}
};
alignas(64) uint8_t pext_table[16][16] = {
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1},
{0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2},
{0, 1, 1, 3, 0, 1, 1, 3, 0, 1, 1, 3, 0, 1, 1, 3},
{0, 0, 0, 0, 1, 2, 2, 4, 0, 0, 0, 0, 1, 2, 2, 4},
{0, 1, 0, 1, 1, 3, 2, 5, 0, 1, 0, 1, 1, 3, 2, 5},
{0, 0, 1, 2, 1, 2, 3, 6, 0, 0, 1, 2, 1, 2, 3, 6},
{0, 1, 1, 3, 1, 3, 3, 7, 0, 1, 1, 3, 1, 3, 3, 7},
{0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2, 4, 2, 4, 4, 8},
{0, 1, 0, 1, 0, 1, 0, 1, 1, 3, 2, 5, 2, 5, 4, 9},
{0, 0, 1, 2, 0, 0, 1, 2, 1, 2, 3, 6, 2, 4, 5, 10},
{0, 1, 1, 3, 0, 1, 1, 3, 1, 3, 3, 7, 2, 5, 5, 11},
{0, 0, 0, 0, 1, 2, 2, 4, 1, 2, 2, 4, 3, 6, 6, 12},
{0, 1, 0, 1, 1, 3, 2, 5, 1, 3, 2, 5, 3, 7, 6, 13},
{0, 0, 1, 2, 1, 2, 3, 6, 1, 2, 3, 6, 3, 6, 7, 14},
{0, 1, 1, 3, 1, 3, 3, 7, 1, 3, 3, 7, 3, 7, 7, 15}
};
alignas(64) uint8_t popcount_table[16] = { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };


inline uint64_t pdep_plain(uint64_t a, uint64_t mask) noexcept {

	uint64_t answer = 0;

	for (uint64_t m = 0; mask; mask >>= 4, m += 4) {
		const uint64_t b = mask % 16;
		answer += ((uint64_t)pdep_table[a % 16][b]) << m;

		a >>= popcount_table[b];
	}

	return answer;
}

inline uint64_t pdep64(uint64_t x, uint64_t mask) noexcept {

#ifdef ONLY_BASIC_INSTRUCTION
	return pdep_plain(x, mask);
#else
	return _pdep_u64(x, mask);
#endif

}

inline uint64_t pext_plain(uint64_t a, uint64_t mask) noexcept {

	uint64_t answer = 0;

	for (uint64_t m = 0; mask; mask >>= 4, a >>= 4) {
		const uint64_t b = mask % 16;
		answer += ((uint64_t)pext_table[a % 16][b]) << m;
		m += popcount_table[b];
	}

	return answer;
}

inline uint64_t pext64(uint64_t x, uint64_t mask) noexcept {

#ifdef ONLY_BASIC_INSTRUCTION
	return pext_plain(x, mask);
#else
	return _pext_u64(x, mask);
#endif

}

/* Original C code from: http://xoroshiro.di.unimi.it/xoroshiro128plus.c */
struct Xoroshiro128plus {

private:

	uint64_t splitmix64(uint64_t &x) noexcept {
		x += 0x9E37'79B9'7F4A'7C15ULL;
		uint64_t z = x;
		z = (z ^ (z >> 30)) * 0xBF58'476D'1CE4'E5B9ULL;
		z = (z ^ (z >> 27)) * 0x94D0'49BB'1331'11EBULL;
		return z ^ (z >> 31);
	}

	uint64_t rotl(const uint64_t x, const int k) noexcept {
		return (x << k) | (x >> (64 - k));
	}

public:

	uint64_t state[2];

	uint64_t rand() noexcept {
		const uint64_t s0 = state[0];
		uint64_t s1 = state[1];
		const uint64_t result = s0 + s1;
		s1 ^= s0;
		state[0] = rotl(s0, 55) ^ s1 ^ (s1 << 14); // a, b
		state[1] = rotl(s1, 36); // c
		return result;
	}

	Xoroshiro128plus() noexcept {
		uint64_t seed = 0;
		state[0] = splitmix64(seed);
		state[1] = splitmix64(seed);
	}

	Xoroshiro128plus(uint64_t seed) noexcept {
		state[0] = splitmix64(seed);
		state[1] = splitmix64(seed);
	}

	Xoroshiro128plus(const uint64_t seed1, const uint64_t seed2) noexcept {
		state[0] = seed1;
		state[1] = seed2;
	}
};

//ハッカーのたのしみ7章124ページ以降で論じられているcompress関数である。x64ではpextという名前で実装されている。
uint64_t compress(uint64_t x, uint64_t mask) noexcept {
	return pext64(x, mask);
}

//ハッカーのたのしみ7章132ページのSAG関数（の亜種）である。pos_m==popcnt(~m)を仮定する。
uint64_t sag(uint64_t x, uint64_t m, uint8_t pos_m) noexcept {
	return (compress(x, m) << pos_m) | compress(x, ~m);
}

uint64_t sag(uint64_t x, uint64_t m) noexcept {
	return (compress(x, m) << (popcount64(~m))) | compress(x, ~m);
}

//ハッカーのたのしみ7章107ページで紹介されている。
uint64_t bitwise_reverse(uint64_t x) noexcept {
	x = ((x & 0x5555'5555'5555'5555ULL) << 1) | ((x >> 1) & 0x5555'5555'5555'5555ULL);
	x = ((x & 0x3333'3333'3333'3333ULL) << 2) | ((x >> 2) & 0x3333'3333'3333'3333ULL);
	x = ((x & 0x0F0F'0F0F'0F0F'0F0FULL) << 4) | ((x >> 4) & 0x0F0F'0F0F'0F0F'0F0FULL);
	x = ((x & 0x00FF'00FF'00FF'00FFULL) << 8) | ((x >> 8) & 0x00FF'00FF'00FF'00FFULL);
	return (x << 48) | ((x << 16) & 0x0000'FFFF'0000'0000ULL) | ((x >> 16) & 0x0000'0000'FFFF'0000ULL) | (x >> 48);
}

//xとyを2進数表示したときに、xが"pattern"になっている箇所で、かつその0の箇所でyのビットが立っているような箇所をすべて探して、
//その最下位ビットが立っている値を返す。ただし、xとyで同じ箇所にビットが立っていないことを仮定する。
uint64_t find_bitpattern_naive(uint64_t x, uint64_t y, uint64_t pattern, uint64_t length) {

	assert((x & y) == 0);
	assert(length <= 64);
	if (length < 64) {
		assert((pattern & ((1ULL << length) - 1ULL)) == pattern);
	}

	const uint64_t neg_pattern = (~pattern) & ((1ULL << length) - 1ULL);

	uint64_t answer = 0;
	for (int i = 0; i < 64; ++i, x >>= 1, y >>= 1) {
		if ((x & pattern) == pattern && (y & neg_pattern) == neg_pattern)answer |= 1ULL << i;
	}
	return answer;
}

//xを2進数表示したときに、1が4個連続している箇所を全て探して、その最下位ビットが立っている値を返す。
//この関数はハッカーのたのしみ106ページ図6.5を参考にして書かれた。
uint64_t find_4_or_more_consecutive_bits(uint64_t x) noexcept {
	x = x & (x >> 2);
	x = x & (x >> 1);
	return x;
}

//xを2進数表示したときに、1が3個連続している箇所を全て探して、その最下位ビットが立っている値を返す。
//この関数はハッカーのたのしみ106ページ図6.5を参考にして書かれた。
uint64_t find_3_or_more_consecutive_bits(uint64_t x) noexcept {
	return x & (x >> 1) & (x >> 2);
}

//xとyを2進数表示したときに、xが"1011"になっている箇所で、かつその0の箇所でyのビットが立っているような箇所をすべて探して、
//その最下位ビットが立っている値を返す。ただし、xとyで同じ箇所にビットが立っていないことを仮定する。
uint64_t find_0b1011(const uint64_t x, const uint64_t y) noexcept {
	const uint64_t xx = (x >> 1) & x;
	return (x >> 3) & (y >> 2) & xx;
}

//xとyを2進数表示したときに、xが"1101"になっている箇所で、かつその0の箇所でyのビットが立っているような箇所をすべて探して、
//その最下位ビットが立っている値を返す。ただし、xとyで同じ箇所にビットが立っていないことを仮定する。
uint64_t find_0b1101(const uint64_t x, const uint64_t y) noexcept {
	const uint64_t xx = (x >> 1) & x;
	return (xx >> 2) & (y >> 1) & x;
}

//xとyを2進数表示したときに、xが"1100"になっている箇所で、かつその0の箇所でyのビットが立っているような箇所をすべて探して、
//その最下位ビットが立っている値を返す。ただし、xとyで同じ箇所にビットが立っていないことを仮定する。
uint64_t find_0b1100(const uint64_t x, const uint64_t y) noexcept {
	const uint64_t xx = (x >> 1) & x;
	const uint64_t yy = (y >> 1) & y;
	return (xx >> 2) & yy;
}

//xとyを2進数表示したときに、xが"0011"になっている箇所で、かつその0の箇所でyのビットが立っているような箇所をすべて探して、
//その最下位ビットが立っている値を返す。ただし、xとyで同じ箇所にビットが立っていないことを仮定する。
uint64_t find_0b0011(const uint64_t x, const uint64_t y) noexcept {
	const uint64_t xx = (x >> 1) & x;
	const uint64_t yy = (y >> 1) & y;
	return (yy >> 2) & xx;
}

//xとyを2進数表示したときに、xが"1001"になっている箇所で、かつその0の箇所でyのビットが立っているような箇所をすべて探して、
//その最下位ビットが立っている値を返す。ただし、xとyで同じ箇所にビットが立っていないことを仮定する。
uint64_t find_0b1001(const uint64_t x, const uint64_t y) noexcept {
	const uint64_t yy = (y >> 1) & y;
	return (x >> 3) & (yy >> 1) & x;
}

//xとyを2進数表示したときに、xが"0101"になっている箇所で、かつその0の箇所でyのビットが立っているような箇所をすべて探して、
//その最下位ビットが立っている値を返す。ただし、xとyで同じ箇所にビットが立っていないことを仮定する。
uint64_t find_0b0101(const uint64_t x, const uint64_t y) noexcept {
	const uint64_t x_x = (x >> 2) & x;
	const uint64_t y_y = (y >> 2) & y;
	return (y_y >> 1) & x_x;
}

//xとyを2進数表示したときに、xが"1010"になっている箇所で、かつその0の箇所でyのビットが立っているような箇所をすべて探して、
//その最下位ビットが立っている値を返す。ただし、xとyで同じ箇所にビットが立っていないことを仮定する。
uint64_t find_0b1010(const uint64_t x, const uint64_t y) noexcept {
	const uint64_t x_x = (x >> 2) & x;
	const uint64_t y_y = (y >> 2) & y;
	return (x_x >> 1) & y_y;
}

//xとyを2進数表示したときに、xが"101"になっている箇所で、かつその0の箇所でyのビットが立っているような箇所をすべて探して、
//その最下位ビットが立っている値を返す。ただし、xとyで同じ箇所にビットが立っていないことを仮定する。
uint64_t find_0b101(const uint64_t x, const uint64_t y) noexcept {
	return (x >> 2) & (y >> 1) & x;
}

//xとyを2進数表示したときに、xが"011"になっている箇所で、かつその0の箇所でyのビットが立っているような箇所をすべて探して、
//その最下位ビットが立っている値を返す。ただし、xとyで同じ箇所にビットが立っていないことを仮定する。
uint64_t find_0b011(const uint64_t x, const uint64_t y) noexcept {
	const uint64_t xx = (x >> 1) & x;
	return (y >> 2) & xx;
}

//xとyを2進数表示したときに、xが"110"になっている箇所で、かつその0の箇所でyのビットが立っているような箇所をすべて探して、
//その最下位ビットが立っている値を返す。ただし、xとyで同じ箇所にビットが立っていないことを仮定する。
uint64_t find_0b110(const uint64_t x, const uint64_t y) noexcept {
	const uint64_t xx = (x >> 1) & x;
	return (xx >> 1) & y;
}

constexpr uint64_t BB_ALL = 0x1FFF'FFFF'FFFF'FFFFULL;

const uint8_t length_table[9] = { 5,6,7,8,9,8,7,6,5 };
const uint8_t start_table[9] = { 0,5,11,18,26,35,43,50,56 };

//uint64_t型変数がYavalathのビットボードだとする。すなわち、
/*
	  1 2 3 4 5
	a . . . . . 6
   b . . . . . . 7
  c . . . . . . . 8
 d . . . . . . . . 9
e . . . . . . . . .
 f . . . . . . . . 9
  g . . . . . . . 8
   h . . . . . . 7
	i . . . . . 6
	  1 2 3 4 5
*/
//という盤面（例えば左上はa1、右下はb5などと呼ぶ）だとして、変数が
//変数の最下位ビットが[a1]で、その次のビットが[a2]で、0b100000のビットが[b1]で、2^60のビットが[i5]だとする。

//文字列をビットボードの座標に変換する。
uint8_t str2pos_bitboard(const std::string pos) {

	const int a1 = pos[0] - 'a';
	if (a1 < 0 || 9 <= a1) {
		throw std::out_of_range("str2pos_bitboard: " + pos);
	}
	const int a2 = pos[1] - '1';
	if (a2 < 0 || length_table[a1] <= a2) {
		throw std::out_of_range("str2pos_bitboard: " + pos);
	}

	return start_table[a1] + a2;
}

//ビットボードの座標を文字列で指定したとして、そこの値を返す。
bool get1pos_bitboard(const uint64_t x, const std::string pos) {
	return (x & (1ULL << str2pos_bitboard(pos))) != 0;
}

//ビットボードの座標を文字列で指定したとして、その座標のビットだけが立っているビットボードを返す。
uint64_t set1pos_bitboard(const std::string pos) {
	return 1ULL << str2pos_bitboard(pos);
}

//ビットボードを時計回りに60度回転した値を求めて返す。
uint64_t rotate_bitboard_060degree_clockwise_naive(const uint64_t x) {

	assert((x & BB_ALL) == x);

	uint64_t answer = 0;

	if (get1pos_bitboard(x, "e1"))answer |= set1pos_bitboard("a1");
	if (get1pos_bitboard(x, "d1"))answer |= set1pos_bitboard("a2");
	if (get1pos_bitboard(x, "c1"))answer |= set1pos_bitboard("a3");
	if (get1pos_bitboard(x, "b1"))answer |= set1pos_bitboard("a4");
	if (get1pos_bitboard(x, "a1"))answer |= set1pos_bitboard("a5");

	if (get1pos_bitboard(x, "f1"))answer |= set1pos_bitboard("b1");
	if (get1pos_bitboard(x, "e2"))answer |= set1pos_bitboard("b2");
	if (get1pos_bitboard(x, "d2"))answer |= set1pos_bitboard("b3");
	if (get1pos_bitboard(x, "c2"))answer |= set1pos_bitboard("b4");
	if (get1pos_bitboard(x, "b2"))answer |= set1pos_bitboard("b5");
	if (get1pos_bitboard(x, "a2"))answer |= set1pos_bitboard("b6");

	if (get1pos_bitboard(x, "g1"))answer |= set1pos_bitboard("c1");
	if (get1pos_bitboard(x, "f2"))answer |= set1pos_bitboard("c2");
	if (get1pos_bitboard(x, "e3"))answer |= set1pos_bitboard("c3");
	if (get1pos_bitboard(x, "d3"))answer |= set1pos_bitboard("c4");
	if (get1pos_bitboard(x, "c3"))answer |= set1pos_bitboard("c5");
	if (get1pos_bitboard(x, "b3"))answer |= set1pos_bitboard("c6");
	if (get1pos_bitboard(x, "a3"))answer |= set1pos_bitboard("c7");

	if (get1pos_bitboard(x, "h1"))answer |= set1pos_bitboard("d1");
	if (get1pos_bitboard(x, "g2"))answer |= set1pos_bitboard("d2");
	if (get1pos_bitboard(x, "f3"))answer |= set1pos_bitboard("d3");
	if (get1pos_bitboard(x, "e4"))answer |= set1pos_bitboard("d4");
	if (get1pos_bitboard(x, "d4"))answer |= set1pos_bitboard("d5");
	if (get1pos_bitboard(x, "c4"))answer |= set1pos_bitboard("d6");
	if (get1pos_bitboard(x, "b4"))answer |= set1pos_bitboard("d7");
	if (get1pos_bitboard(x, "a4"))answer |= set1pos_bitboard("d8");

	if (get1pos_bitboard(x, "i1"))answer |= set1pos_bitboard("e1");
	if (get1pos_bitboard(x, "h2"))answer |= set1pos_bitboard("e2");
	if (get1pos_bitboard(x, "g3"))answer |= set1pos_bitboard("e3");
	if (get1pos_bitboard(x, "f4"))answer |= set1pos_bitboard("e4");
	if (get1pos_bitboard(x, "e5"))answer |= set1pos_bitboard("e5");
	if (get1pos_bitboard(x, "d5"))answer |= set1pos_bitboard("e6");
	if (get1pos_bitboard(x, "c5"))answer |= set1pos_bitboard("e7");
	if (get1pos_bitboard(x, "b5"))answer |= set1pos_bitboard("e8");
	if (get1pos_bitboard(x, "a5"))answer |= set1pos_bitboard("e9");

	if (get1pos_bitboard(x, "i2"))answer |= set1pos_bitboard("f1");
	if (get1pos_bitboard(x, "h3"))answer |= set1pos_bitboard("f2");
	if (get1pos_bitboard(x, "g4"))answer |= set1pos_bitboard("f3");
	if (get1pos_bitboard(x, "f5"))answer |= set1pos_bitboard("f4");
	if (get1pos_bitboard(x, "e6"))answer |= set1pos_bitboard("f5");
	if (get1pos_bitboard(x, "d6"))answer |= set1pos_bitboard("f6");
	if (get1pos_bitboard(x, "c6"))answer |= set1pos_bitboard("f7");
	if (get1pos_bitboard(x, "b6"))answer |= set1pos_bitboard("f8");

	if (get1pos_bitboard(x, "i3"))answer |= set1pos_bitboard("g1");
	if (get1pos_bitboard(x, "h4"))answer |= set1pos_bitboard("g2");
	if (get1pos_bitboard(x, "g5"))answer |= set1pos_bitboard("g3");
	if (get1pos_bitboard(x, "f6"))answer |= set1pos_bitboard("g4");
	if (get1pos_bitboard(x, "e7"))answer |= set1pos_bitboard("g5");
	if (get1pos_bitboard(x, "d7"))answer |= set1pos_bitboard("g6");
	if (get1pos_bitboard(x, "c7"))answer |= set1pos_bitboard("g7");

	if (get1pos_bitboard(x, "i4"))answer |= set1pos_bitboard("h1");
	if (get1pos_bitboard(x, "h5"))answer |= set1pos_bitboard("h2");
	if (get1pos_bitboard(x, "g6"))answer |= set1pos_bitboard("h3");
	if (get1pos_bitboard(x, "f7"))answer |= set1pos_bitboard("h4");
	if (get1pos_bitboard(x, "e8"))answer |= set1pos_bitboard("h5");
	if (get1pos_bitboard(x, "d8"))answer |= set1pos_bitboard("h6");

	if (get1pos_bitboard(x, "i5"))answer |= set1pos_bitboard("i1");
	if (get1pos_bitboard(x, "h6"))answer |= set1pos_bitboard("i2");
	if (get1pos_bitboard(x, "g7"))answer |= set1pos_bitboard("i3");
	if (get1pos_bitboard(x, "f8"))answer |= set1pos_bitboard("i4");
	if (get1pos_bitboard(x, "e9"))answer |= set1pos_bitboard("i5");

	return answer;
}

//ビットボードを時計回りに120度回転した値を求めて返す。
uint64_t rotate_bitboard_120degree_clockwise_naive(const uint64_t x) {

	assert((x & BB_ALL) == x);

	return rotate_bitboard_060degree_clockwise_naive(rotate_bitboard_060degree_clockwise_naive(x));
}

//ビットボードを時計回りに180度回転した値を求めて返す。
uint64_t rotate_bitboard_180degree_clockwise_naive(const uint64_t x) {

	assert((x & BB_ALL) == x);

	return rotate_bitboard_060degree_clockwise_naive(rotate_bitboard_120degree_clockwise_naive(x));
}

//ビットボードを時計回りに240度回転した値を求めて返す。
uint64_t rotate_bitboard_240degree_clockwise_naive(const uint64_t x) {

	assert((x & BB_ALL) == x);

	return rotate_bitboard_060degree_clockwise_naive(rotate_bitboard_180degree_clockwise_naive(x));
}

//ビットボードを時計回りに300度回転した値を求めて返す。
uint64_t rotate_bitboard_300degree_clockwise_naive(const uint64_t x) {

	assert((x & BB_ALL) == x);

	return rotate_bitboard_060degree_clockwise_naive(rotate_bitboard_240degree_clockwise_naive(x));
}

//ビットボードを横向きに反転した値を求めて返す。
uint64_t mirror_bitboard_naive(const uint64_t x) {

	assert((x & BB_ALL) == x);

	uint64_t answer = 0;

	for (int i = 0; i < 9; ++i) {
		for (int j = 0; j < length_table[i]; ++j) {
			const uint64_t b1 = 1ULL << (start_table[i] + j);
			const uint64_t b2 = 1ULL << (start_table[i] + length_table[i] - 1 - j);
			if (x & b1)answer |= b2;
		}
	}

	return answer;
}

uint8_t rot060_pos[61] = {}, rot120_pos[61] = {}, rot240_pos[61] = {}, rot300_pos[61] = {}, mirror_pos[61] = {};

void init_rot_pos_arrays(uint8_t __rot060_pos[61], uint8_t __rot120_pos[61], uint8_t __rot240_pos[61], uint8_t __rot300_pos[61], uint8_t __mirror_pos[61]) {

	for (int i = 0; i < 61; ++i) {
		const uint64_t src_1bit = 1ULL << i;

		const uint64_t dst_1bit_060degree = rotate_bitboard_060degree_clockwise_naive(src_1bit);
		unsigned long dst_pos_060degree = 0;
		bitscan_forward64(&dst_pos_060degree, dst_1bit_060degree);
		assert(dst_pos_060degree < 61);
		__rot060_pos[i] = uint8_t(dst_pos_060degree);

		const uint64_t dst_1bit_120degree = rotate_bitboard_120degree_clockwise_naive(src_1bit);
		unsigned long dst_pos_120degree = 0;
		bitscan_forward64(&dst_pos_120degree, dst_1bit_120degree);
		assert(dst_pos_120degree < 61);
		__rot120_pos[i] = uint8_t(dst_pos_120degree);

		const uint64_t dst_1bit_240degree = rotate_bitboard_240degree_clockwise_naive(src_1bit);
		unsigned long dst_pos_240degree = 0;
		bitscan_forward64(&dst_pos_240degree, dst_1bit_240degree);
		assert(dst_pos_240degree < 61);
		__rot240_pos[i] = uint8_t(dst_pos_240degree);

		const uint64_t dst_1bit_300degree = rotate_bitboard_300degree_clockwise_naive(src_1bit);
		unsigned long dst_pos_300degree = 0;
		bitscan_forward64(&dst_pos_300degree, dst_1bit_300degree);
		assert(dst_pos_300degree < 61);
		__rot300_pos[i] = uint8_t(dst_pos_300degree);

		const uint64_t dst_1bit_mirror = mirror_bitboard_naive(src_1bit);
		unsigned long dst_pos_mirror = 0;
		bitscan_forward64(&dst_pos_mirror, dst_1bit_mirror);
		assert(dst_pos_mirror < 61);
		__mirror_pos[i] = uint8_t(dst_pos_mirror);
	}
}

//ハッカーのたのしみ132ページのp配列。
uint64_t p_rot060[6] = {}, p_rot120[6] = {}, p_rot240[6] = {}, p_rot300[6] = {}, p_mirror[6] = {};
uint8_t p_rot060_pop[6] = {}, p_rot120_pop[6] = {}, p_rot240_pop[6] = {}, p_rot300_pop[6] = {}, p_mirror_pop[6] = {};






void init_p_arrays___(
	uint64_t __p_rot060[6], uint64_t __p_rot120[6], uint64_t __p_rot240[6], uint64_t __p_rot300[6], uint64_t __p_mirror[6],
	uint8_t __p_rot060_pop[6], uint8_t __p_rot120_pop[6], uint8_t __p_rot240_pop[6], uint8_t __p_rot300_pop[6], uint8_t __p_mirror_pop[6]) {

	//ハッカーのたのしみ133ページの事前計算を行う。
	const auto sagfunc = [](uint64_t p[6]) {
		p[1] = sag(p[1], p[0]);
		p[2] = sag(sag(p[2], p[0]), p[1]);
		p[3] = sag(sag(sag(p[3], p[0]), p[1]), p[2]);
		p[4] = sag(sag(sag(sag(p[4], p[0]), p[1]), p[2]), p[3]);
		p[5] = sag(sag(sag(sag(sag(p[5], p[0]), p[1]), p[2]), p[3]), p[4]);
	};

	const auto func = [sagfunc](uint64_t __p[6], uint8_t __pos[61], uint8_t __p_pop[6]) {
		for (int i = 0; i < 6; ++i)__p[i] = 0;
		for (int i = 0; i < 61; ++i) {
			for (uint64_t b = __pos[i], j = 0; b; ++j, b >>= 1) {
				if (b & 1ULL)__p[j] |= 1ULL << i;
			}
		}
		for (uint64_t i = 61; i < 64; ++i) {
			for (uint64_t b = i, j = 0; b; ++j, b >>= 1) {
				if (b & 1ULL)__p[j] |= 1ULL << i;
			}
		}
		sagfunc(__p);
		for (int i = 0; i < 6; ++i)__p_pop[i] = popcount64(__p[i]);
	};

	uint8_t __rot060_pos[61] = {}, __rot120_pos[61] = {}, __rot240_pos[61] = {}, __rot300_pos[61] = {}, __mirror_pos[61] = {};
	init_rot_pos_arrays(__rot060_pos, __rot120_pos, __rot240_pos, __rot300_pos, __mirror_pos);

	func(__p_rot060, __rot060_pos, __p_rot060_pop);
	func(__p_rot120, __rot120_pos, __p_rot120_pop);
	func(__p_rot240, __rot240_pos, __p_rot240_pop);
	func(__p_rot300, __rot300_pos, __p_rot300_pop);
	func(__p_mirror, __mirror_pos, __p_mirror_pop);
}

void init_p_arrays(
	uint64_t __p_rot060[6], uint64_t __p_rot120[6], uint64_t __p_rot240[6], uint64_t __p_rot300[6], uint64_t __p_mirror[6],
	uint8_t __p_rot060_pop[6], uint8_t __p_rot120_pop[6], uint8_t __p_rot240_pop[6], uint8_t __p_rot300_pop[6], uint8_t __p_mirror_pop[6]) {

	//ハッカーのたのしみ133ページの事前計算を行う。
	const auto sagfunc = [](uint64_t p[6]) {
		p[1] = sag(p[1], p[0]);
		p[2] = sag(sag(p[2], p[0]), p[1]);
		p[3] = sag(sag(sag(p[3], p[0]), p[1]), p[2]);
		p[4] = sag(sag(sag(sag(p[4], p[0]), p[1]), p[2]), p[3]);
		p[5] = sag(sag(sag(sag(sag(p[5], p[0]), p[1]), p[2]), p[3]), p[4]);
	};

	const auto func = [sagfunc](uint64_t __p[6], uint8_t __pos[61], uint8_t __p_pop[6]) {
		for (int i = 0; i < 6; ++i)__p[i] = 0;

		uint8_t pos64[64];
		for (int i = 0; i < 61; ++i)pos64[i] = __pos[i];
		for (int i = 61; i < 64; ++i)pos64[i] = i;

		uint8_t pos_stable64[64];

		int rep = 0;
		for (int count = 0; count < 64; ++rep) {
			for (int i = 0; i < 64; ++i) {
				if (pos64[i] == count) {
					pos_stable64[i] = rep;
					count++;
				}
			}
		}

		for (int i = 0; i < 64; ++i) {
			for (uint64_t b = pos_stable64[i], j = 0; b; ++j, b >>= 1) {
				if (b & 1ULL) {
					__p[j] |= 1ULL << i;
				}
			}
		}

		sagfunc(__p);

		for (int i = 0; i < 6; ++i)__p_pop[i] = 64 - popcount64(__p[i]);
	};

	uint8_t __rot060_pos[61] = {}, __rot120_pos[61] = {}, __rot240_pos[61] = {}, __rot300_pos[61] = {}, __mirror_pos[61] = {};
	init_rot_pos_arrays(__rot060_pos, __rot120_pos, __rot240_pos, __rot300_pos, __mirror_pos);

	func(__p_rot060, __rot060_pos, __p_rot060_pop);
	func(__p_rot120, __rot120_pos, __p_rot120_pop);
	func(__p_rot240, __rot240_pos, __p_rot240_pop);
	func(__p_rot300, __rot300_pos, __p_rot300_pop);
	func(__p_mirror, __mirror_pos, __p_mirror_pop);
}







//ビットボードを時計回りに60度回転した値を求めて返す。
uint64_t rotate_bitboard_060degree_clockwise(uint64_t x) noexcept {

	x = sag(x, p_rot060[0], p_rot060_pop[0]);
	x = sag(x, p_rot060[1], p_rot060_pop[1]);
	x = sag(x, p_rot060[2], p_rot060_pop[2]);
	x = sag(x, p_rot060[3], p_rot060_pop[3]);
	x = sag(x, p_rot060[4], p_rot060_pop[4]);
	x = sag(x, p_rot060[5], p_rot060_pop[5]);
	return x;
}

//ビットボードを時計回りに120度回転した値を求めて返す。
uint64_t rotate_bitboard_120degree_clockwise(uint64_t x) noexcept {

	x = sag(x, p_rot120[0], p_rot120_pop[0]);
	x = sag(x, p_rot120[1], p_rot120_pop[1]);
	x = sag(x, p_rot120[2], p_rot120_pop[2]);
	x = sag(x, p_rot120[3], p_rot120_pop[3]);
	x = sag(x, p_rot120[4], p_rot120_pop[4]);
	x = sag(x, p_rot120[5], p_rot120_pop[5]);
	return x;
}

//ビットボードを時計回りに180度回転した値を求めて返す。
uint64_t rotate_bitboard_180degree_clockwise(uint64_t x) noexcept {

	return bitwise_reverse(x) >> 3;
}

//ビットボードを時計回りに240度回転した値を求めて返す。
uint64_t rotate_bitboard_240degree_clockwise(uint64_t x) noexcept {

	x = sag(x, p_rot240[0], p_rot240_pop[0]);
	x = sag(x, p_rot240[1], p_rot240_pop[1]);
	x = sag(x, p_rot240[2], p_rot240_pop[2]);
	x = sag(x, p_rot240[3], p_rot240_pop[3]);
	//x = sag(x, p_rot240[4], p_rot240_pop[4]);
	//x = sag(x, p_rot240[5], p_rot240_pop[5]);
	return x;
}

//ビットボードを時計回りに300度回転した値を求めて返す。
uint64_t rotate_bitboard_300degree_clockwise(uint64_t x) noexcept {

	x = sag(x, p_rot300[0], p_rot300_pop[0]);
	x = sag(x, p_rot300[1], p_rot300_pop[1]);
	x = sag(x, p_rot300[2], p_rot300_pop[2]);
	x = sag(x, p_rot300[3], p_rot300_pop[3]);
	//x = sag(x, p_rot300[4], p_rot300_pop[4]);
	//x = sag(x, p_rot300[5], p_rot300_pop[5]);
	return x;
}

//ビットボードを横向きに反転した値を求めて返す。
uint64_t mirror_bitboard(uint64_t x) noexcept {

	x = sag(x, p_mirror[0], p_mirror_pop[0]);
	x = sag(x, p_mirror[1], p_mirror_pop[1]);
	x = sag(x, p_mirror[2], p_mirror_pop[2]);
	x = sag(x, p_mirror[3], p_mirror_pop[3]);
	x = sag(x, p_mirror[4], p_mirror_pop[4]);
	x = sag(x, p_mirror[5], p_mirror_pop[5]);
	return x;
}

constexpr uint64_t BB_EDGE_2 = 0b000'11111'111111'1100011'11000011'110000011'11000011'1100011'111111'11111ULL;
constexpr uint64_t BB_EDGE_3 = 0b000'11111'111111'1111111'11100111'111000111'11100111'1111111'111111'11111ULL;

constexpr uint64_t BB_MASK_4 = 0b000'00011'000111'0001111'00011111'000111111'00011111'0001111'000111'00011ULL;
constexpr uint64_t BB_MASK_3 = 0b000'00111'001111'0011111'00111111'001111111'00111111'0011111'001111'00111ULL;
constexpr uint64_t BB_MASK_2 = 0b000'01111'011111'0111111'01111111'011111111'01111111'0111111'011111'01111ULL;

constexpr uint64_t BB_MASK_L = 0b000'00001'000001'0000001'00000001'000000001'00000001'0000001'000001'00001ULL;
constexpr uint64_t BB_MASK_R = 0b000'10000'100000'1000000'10000000'100000000'10000000'1000000'100000'10000ULL;

constexpr uint64_t BB_MASK_R2 = 0b000'11110'111110'1111110'11111110'111111110'11111110'1111110'111110'11110ULL;
constexpr uint64_t BB_MASK_R3 = 0b000'11100'111100'1111100'11111100'111111100'11111100'1111100'111100'11100ULL;

struct Bitboard {

public:
	uint64_t b000, b060, b120;

	Bitboard() noexcept {
		b000 = 0;
		b060 = 0;
		b120 = 0;
	}

	Bitboard(const uint64_t _b000) noexcept {
		b000 = _b000;
		b060 = rotate_bitboard_060degree_clockwise(_b000);
		b120 = rotate_bitboard_120degree_clockwise(_b000);
	}

	bool validate() const noexcept {
		if ((b000 & BB_ALL) != b000)return false;
		if (rotate_bitboard_060degree_clockwise(b000) != b060)return false;
		if (rotate_bitboard_120degree_clockwise(b000) != b120)return false;
		return true;
	}

	bool test(const uint8_t pos) const noexcept {
		assert(pos < 61);
		assert(validate());
		return (b000 & (1ULL << pos)) != 0;
	}

	void do_move(const uint8_t pos) noexcept {
		assert(pos < 61);
		assert(validate());
		assert((b000 & (1ULL << pos)) == 0);
		b000 |= 1ULL << pos;
		b060 |= 1ULL << rot060_pos[pos];
		b120 |= 1ULL << rot120_pos[pos];
	}

	void undo_move(const uint8_t pos) noexcept {
		assert(pos < 61);
		assert(validate());
		assert((b000 & (1ULL << pos)) != 0);
		b000 ^= 1ULL << pos;
		b060 ^= 1ULL << rot060_pos[pos];
		b120 ^= 1ULL << rot120_pos[pos];
	}
};

void print_board(const uint64_t player, const uint64_t opponent) {
	/*
		  1 2 3 4 5
		a . . . . . 6
	   b . . . . . . 7
	  c . . . . . . . 8
	 d . . . . . . . . 9
	e . . . . . . . . .
	 f . . . . . . . . 9
	  g . . . . . . . 8
	   h . . . . . . 7
		i . . . . . 6
		  1 2 3 4 5
	*/

	assert((player & opponent) == 0);
	assert(player <= BB_ALL);
	assert(opponent <= BB_ALL);

	const int len_space[9] = { 4,3,2,1,0,1,2,3,4 };
	const std::string L = "abcdefghi";
	const std::string R = "6789 9876";

	std::cout << "      1 2 3 4 5" << std::endl;
	for (int i = 0; i < 9; ++i) {
		for (int j = 0; j < len_space[i]; ++j)std::cout << " ";
		std::cout << L[i] << " ";
		for (int j = 0; j < length_table[i]; ++j) {
			const uint64_t bit = 1ULL << (start_table[i] + j);
			if (player & bit) {
				std::cout << "X ";
			}
			else if (opponent & bit) {
				std::cout << "O ";
			}
			else {
				std::cout << ". ";
			}
		}
		std::cout << R[i] << std::endl;
	}
	std::cout << "      1 2 3 4 5" << std::endl;
}
void print_board(const std::array<int, 62>&kifu) {

	uint64_t player = 0, opponent = 0;

	for (int i = 0; kifu[i] != -1; ++i) {
		assert(0 <= kifu[i] && kifu[i] < 61);
		if (i % 2 == 0) {
			assert((player & (1ULL << kifu[i])) == 0);
			player |= 1ULL << kifu[i];
		}
		else {
			assert((opponent & (1ULL << kifu[i])) == 0);
			opponent |= 1ULL << kifu[i];
		}
	}

	print_board(player, opponent);
}

//勝利か敗北かを判定する。
enum Immediate {
	UNKNOWN = 1,
	WIN = 1000, //勝利(盤上の駒が4つ以上連続している箇所が存在する)
	LOSE = -1000, //敗北(勝利条件を満たしておらず、かつ盤上の駒が3つ以上連続している箇所が存在する)
	DRAW = 0//引き分け(盤上に空きマスが存在しない)
};
int is_immediate_end(const Bitboard &player) noexcept {
	if ((find_4_or_more_consecutive_bits(player.b000) | find_4_or_more_consecutive_bits(player.b060) | find_4_or_more_consecutive_bits(player.b120)) & BB_MASK_4)return Immediate::WIN;
	if ((find_3_or_more_consecutive_bits(player.b000) | find_3_or_more_consecutive_bits(player.b060) | find_3_or_more_consecutive_bits(player.b120)) & BB_MASK_3)return Immediate::LOSE;
	return Immediate::UNKNOWN;
}

bool lessthan_bitboard(const uint64_t a_player, const uint64_t a_opponent, const uint64_t b_player, const uint64_t b_opponent) noexcept {
	//a = a_player * 2^64 + a_opponent, 
	//b = b_player * 2^64 + b_opponent として、a<bならtrueを返す。
	return (a_player < b_player) || ((a_player == b_player && a_opponent < b_opponent));
}
void get_unique_bitboard(const uint64_t player, const uint64_t opponent, uint64_t *result_player, uint64_t *result_opponent) noexcept {

	*result_player = player;
	*result_opponent = opponent;

	//このfor文は全部アンロールしてmirrorの結果とか再利用すれば高速化できる可能性がある。
	for (int i = 1; i < 12; ++i) {
		uint64_t tmp_p = player;
		uint64_t tmp_o = opponent;
		if (i & 1) {
			tmp_p = mirror_bitboard(tmp_p);
			tmp_o = mirror_bitboard(tmp_o);
		}
		if (i & 2) {
			tmp_p = rotate_bitboard_180degree_clockwise(tmp_p);
			tmp_o = rotate_bitboard_180degree_clockwise(tmp_o);
		}
		switch (i / 4) {
		case 0:
			break;
		case 1:
			tmp_p = rotate_bitboard_060degree_clockwise(tmp_p);
			tmp_o = rotate_bitboard_060degree_clockwise(tmp_o);
			break;
		case 2:
			tmp_p = rotate_bitboard_120degree_clockwise(tmp_p);
			tmp_o = rotate_bitboard_120degree_clockwise(tmp_o);
			break;
		default:
			assert(0);
			break;
		}
		if (lessthan_bitboard(tmp_p, tmp_o, *result_player, *result_opponent)) {
			*result_player = tmp_p;
			*result_opponent = tmp_o;
		}
	}
}

//playerがその座標に石を置くと勝ちになるような座標のビットボードを求めて返す。
uint64_t get_winning_pos_bitboard(const Bitboard &player, const Bitboard &empty) noexcept {

	const uint64_t answer000 =
		((find_0b1011(player.b000, empty.b000) & BB_MASK_4) << 2) |
		((find_0b1101(player.b000, empty.b000) & BB_MASK_4) << 1);
	const uint64_t answer060 =
		((find_0b1011(player.b060, empty.b060) & BB_MASK_4) << 2) |
		((find_0b1101(player.b060, empty.b060) & BB_MASK_4) << 1);
	const uint64_t answer120 =
		((find_0b1011(player.b120, empty.b120) & BB_MASK_4) << 2) |
		((find_0b1101(player.b120, empty.b120) & BB_MASK_4) << 1);

	return answer000 | rotate_bitboard_300degree_clockwise(answer060) | rotate_bitboard_240degree_clockwise(answer120);
}

//playerがその座標に石を置くと勝ちになるような座標が存在するならtrue、さもなくばfalseを返す。
bool exist_winning_pos_bitboard(const Bitboard &player, const Bitboard &empty) noexcept {

	const uint64_t answer000 =
		((find_0b1011(player.b000, empty.b000) & BB_MASK_4) << 2) |
		((find_0b1101(player.b000, empty.b000) & BB_MASK_4) << 1);
	const uint64_t answer060 =
		((find_0b1011(player.b060, empty.b060) & BB_MASK_4) << 2) |
		((find_0b1101(player.b060, empty.b060) & BB_MASK_4) << 1);
	const uint64_t answer120 =
		((find_0b1011(player.b120, empty.b120) & BB_MASK_4) << 2) |
		((find_0b1101(player.b120, empty.b120) & BB_MASK_4) << 1);

	return (answer000 | answer060 | answer120) != 0;
}

//playerがその座標に石を置くと3以上連続で石が並ぶような座標のビットボードを求めて返す。
//（負けとイコールではない。4以上連続で石が並ぶなら勝利が優先される。）
uint64_t get_losing_pos_bitboard(const Bitboard &player, const Bitboard &empty) noexcept {

	const uint64_t answer000 =
		((find_0b011(player.b000, empty.b000) & BB_MASK_3) << 2) |
		((find_0b101(player.b000, empty.b000) & BB_MASK_3) << 1) |
		(find_0b110(player.b000, empty.b000) & BB_MASK_3);
	const uint64_t answer060 =
		((find_0b011(player.b060, empty.b060) & BB_MASK_3) << 2) |
		((find_0b101(player.b060, empty.b060) & BB_MASK_3) << 1) |
		(find_0b110(player.b060, empty.b060) & BB_MASK_3);
	const uint64_t answer120 =
		((find_0b011(player.b120, empty.b120) & BB_MASK_3) << 2) |
		((find_0b101(player.b120, empty.b120) & BB_MASK_3) << 1) |
		(find_0b110(player.b120, empty.b120) & BB_MASK_3);

	return answer000 | rotate_bitboard_300degree_clockwise(answer060) | rotate_bitboard_240degree_clockwise(answer120);
}

//playerがその座標に石を置くと王手になるような座標のビットボードを求めて返す。
//（置くと即負けになる座標も含まれる。）
uint64_t get_check_pos_bitboard(const Bitboard &player, const Bitboard &empty) noexcept {

	const uint64_t _0b1001_000 = find_0b1001(player.b000, empty.b000) & BB_MASK_4;
	const uint64_t answer000 =
		(((find_0b0011(player.b000, empty.b000) | find_0b0101(player.b000, empty.b000)) & BB_MASK_4) << 3) |
		(_0b1001_000 << 1) | (_0b1001_000 << 2) |
		((find_0b1100(player.b000, empty.b000) | find_0b1010(player.b000, empty.b000)) & BB_MASK_4);

	const uint64_t _0b1001_060 = find_0b1001(player.b060, empty.b060) & BB_MASK_4;
	const uint64_t answer060 =
		(((find_0b0011(player.b060, empty.b060) | find_0b0101(player.b060, empty.b060)) & BB_MASK_4) << 3) |
		(_0b1001_060 << 1) | (_0b1001_060 << 2) |
		((find_0b1100(player.b060, empty.b060) | find_0b1010(player.b060, empty.b060)) & BB_MASK_4);

	const uint64_t _0b1001_120 = find_0b1001(player.b120, empty.b120) & BB_MASK_4;
	const uint64_t answer120 =
		(((find_0b0011(player.b120, empty.b120) | find_0b0101(player.b120, empty.b120)) & BB_MASK_4) << 3) |
		(_0b1001_120 << 1) | (_0b1001_120 << 2) |
		((find_0b1100(player.b120, empty.b120) | find_0b1010(player.b120, empty.b120)) & BB_MASK_4);

	return answer000 | rotate_bitboard_300degree_clockwise(answer060) | rotate_bitboard_240degree_clockwise(answer120);
}

//lethal_opponent_emptyは、opponentがそこに石を置くと即座に負けになる空白座標のビットボード。
//safe_player_emptyは、playerがそこに石を置いても即座に負けにならない空白座標のビットボード。
//横列にのみ着目するとして、playerが王手をかける手であって、かつopponentがそれを防ぐと即座に3目並んで負けになるような手が存在する⇔trueを返す。
uint64_t can_mate_1ply_1dir(const uint64_t player, const uint64_t opponent, const uint64_t safe_player_empty, const uint64_t lethal_opponent_empty) noexcept {

	const uint64_t bb_p_p = (player >> 1) & player;
	const uint64_t bb_p_x_p = (player >> 2) & player;
	const uint64_t bb_p_x_x_p = (player >> 3) & player;
	const uint64_t bb_s_l_p_p = ((safe_player_empty >> 3) & (lethal_opponent_empty >> 2) & bb_p_p);
	const uint64_t bb_p_p_l_s = ((bb_p_p >> 2) & (lethal_opponent_empty >> 1) & safe_player_empty);
	const uint64_t bb_p_s_l_p = ((safe_player_empty >> 2) & (lethal_opponent_empty >> 1) & bb_p_x_x_p);
	const uint64_t bb_p_l_s_p = ((safe_player_empty >> 1) & (lethal_opponent_empty >> 2) & bb_p_x_x_p);
	const uint64_t bb_p_l_p_s = ((bb_p_x_p >> 1) & (lethal_opponent_empty >> 2) & safe_player_empty);
	const uint64_t bb_s_p_l_p = ((safe_player_empty >> 3) & (lethal_opponent_empty >> 1) & bb_p_x_p);

	return (bb_s_l_p_p | bb_p_p_l_s | bb_p_s_l_p | bb_p_l_s_p | bb_p_l_p_s | bb_s_p_l_p) & BB_MASK_4;
}

//一手詰めかどうかの判定を行う。この局面がplayerの手番だとする。
//すなわち、playerに王手をかける手であって、かつopponentがそれを防ぐと即座に3目並んで負けになるような手が存在する（かつ現在playerに王手がかかっていない）⇔Immediate::WINを返す。
//それに加えて、playerが石を置くと即座に負けになる座標のビットボードを求めて、*player_losing_bitboard_000に代入して返す。
//playerがどこに置いても即座に負けになるか、二重王手がかかっている⇔Immediate::LOSEを返す。さもなくばImmediate::UNKNOWNを返す。
//
//このターンに石を置いて即座に勝利できる場合、この関数は正確な値を返さない。exist_winning_pos_bitboard関数で先にチェックしておくこと。
//
//返り値がImmediate::WINである⇒playerの勝ちである。
//返り値がImmediate::LOSEである⇒playerの負けである。
int search_mate_1ply(const Bitboard &player, const Bitboard &opponent, const Bitboard &empty, const uint64_t opponent_winning, uint64_t &player_losing_bitboard_000) noexcept {

	const uint64_t player_losing_000 =
		((find_0b011(player.b000, empty.b000) & BB_MASK_3) << 2) |
		((find_0b101(player.b000, empty.b000) & BB_MASK_3) << 1) |
		(find_0b110(player.b000, empty.b000) & BB_MASK_3);
	const uint64_t player_losing_060 =
		((find_0b011(player.b060, empty.b060) & BB_MASK_3) << 2) |
		((find_0b101(player.b060, empty.b060) & BB_MASK_3) << 1) |
		(find_0b110(player.b060, empty.b060) & BB_MASK_3);
	const uint64_t player_losing_120 =
		((find_0b011(player.b120, empty.b120) & BB_MASK_3) << 2) |
		((find_0b101(player.b120, empty.b120) & BB_MASK_3) << 1) |
		(find_0b110(player.b120, empty.b120) & BB_MASK_3);
	const uint64_t player_losing = player_losing_000 | rotate_bitboard_300degree_clockwise(player_losing_060) | rotate_bitboard_240degree_clockwise(player_losing_120);

	player_losing_bitboard_000 = player_losing;
	if (player_losing == empty.b000)return Immediate::LOSE;//1

	if (opponent_winning) {
		if ((opponent_winning & player_losing) == opponent_winning)return Immediate::LOSE;//1
		if (opponent_winning & (opponent_winning - 1))return Immediate::LOSE;//2

		//本当は「逆王手の一手詰め」が存在するかもしれないが、そこは調べない。
		//（レアケースだろうし、手を進めればすぐに判明するので、コストに見合わない気がする）
		return Immediate::UNKNOWN;
	}

	const uint64_t opponent_losing_000 =
		((find_0b011(opponent.b000, empty.b000) & BB_MASK_3) << 2) |
		((find_0b101(opponent.b000, empty.b000) & BB_MASK_3) << 1) |
		(find_0b110(opponent.b000, empty.b000) & BB_MASK_3);
	const uint64_t opponent_losing_060 =
		((find_0b011(opponent.b060, empty.b060) & BB_MASK_3) << 2) |
		((find_0b101(opponent.b060, empty.b060) & BB_MASK_3) << 1) |
		(find_0b110(opponent.b060, empty.b060) & BB_MASK_3);
	const uint64_t opponent_losing_120 =
		((find_0b011(opponent.b120, empty.b120) & BB_MASK_3) << 2) |
		((find_0b101(opponent.b120, empty.b120) & BB_MASK_3) << 1) |
		(find_0b110(opponent.b120, empty.b120) & BB_MASK_3);

	if (can_mate_1ply_1dir(player.b000, opponent.b000, empty.b000 & (~player_losing_000), opponent_losing_000))return Immediate::WIN;//2
	if (can_mate_1ply_1dir(player.b060, opponent.b060, empty.b060 & (~player_losing_060), opponent_losing_060))return Immediate::WIN;//2
	if (can_mate_1ply_1dir(player.b120, opponent.b120, empty.b120 & (~player_losing_120), opponent_losing_120))return Immediate::WIN;//2

	return Immediate::UNKNOWN;
}

int search_mate_1ply(const Bitboard &player, const Bitboard &opponent, const Bitboard &empty, const uint64_t opponent_winning, uint64_t &player_losing_bitboard_000, int& bestmove) noexcept {

	const uint64_t player_losing_000 =
		((find_0b011(player.b000, empty.b000) & BB_MASK_3) << 2) |
		((find_0b101(player.b000, empty.b000) & BB_MASK_3) << 1) |
		(find_0b110(player.b000, empty.b000) & BB_MASK_3);
	const uint64_t player_losing_060 =
		((find_0b011(player.b060, empty.b060) & BB_MASK_3) << 2) |
		((find_0b101(player.b060, empty.b060) & BB_MASK_3) << 1) |
		(find_0b110(player.b060, empty.b060) & BB_MASK_3);
	const uint64_t player_losing_120 =
		((find_0b011(player.b120, empty.b120) & BB_MASK_3) << 2) |
		((find_0b101(player.b120, empty.b120) & BB_MASK_3) << 1) |
		(find_0b110(player.b120, empty.b120) & BB_MASK_3);
	const uint64_t player_losing = player_losing_000 | rotate_bitboard_300degree_clockwise(player_losing_060) | rotate_bitboard_240degree_clockwise(player_losing_120);

	player_losing_bitboard_000 = player_losing;

	if (player_losing == empty.b000) {
		unsigned long pos = 0;
		bitscan_forward64(&pos, player_losing);
		bestmove = int(pos);
		return Immediate::LOSE;
	}

	if (opponent_winning) {
		if ((opponent_winning & player_losing) == opponent_winning) {
			unsigned long pos = 0;
			bitscan_forward64(&pos, opponent_winning);
			bestmove = int(pos);
			return Immediate::LOSE;
		}

		if (opponent_winning & (opponent_winning - 1)) {
			const uint64_t bb = opponent_winning & (~player_losing);
			assert(bb);
			unsigned long pos = 0;
			bitscan_forward64(&pos, bb);
			bestmove = int(pos);
			return Immediate::LOSE;
		}

		//本当は「逆王手の一手詰め」が存在するかもしれないが、そこは調べない。
		//（レアケースだろうし、手を進めればすぐに判明するので、コストに見合わない気がする）
		unsigned long pos = 0;
		bitscan_forward64(&pos, opponent_winning);
		bestmove = int(pos);
		return Immediate::UNKNOWN;
	}

	const uint64_t opponent_losing_000 =
		((find_0b011(opponent.b000, empty.b000) & BB_MASK_3) << 2) |
		((find_0b101(opponent.b000, empty.b000) & BB_MASK_3) << 1) |
		(find_0b110(opponent.b000, empty.b000) & BB_MASK_3);
	const uint64_t b000 = can_mate_1ply_1dir(player.b000, opponent.b000, empty.b000 & (~player_losing_000), opponent_losing_000);
	if (b000) {
		unsigned long pos = 0;
		bitscan_forward64(&pos, b000);
		bestmove = int(pos);
		return Immediate::WIN;
	}

	const uint64_t opponent_losing_060 =
		((find_0b011(opponent.b060, empty.b060) & BB_MASK_3) << 2) |
		((find_0b101(opponent.b060, empty.b060) & BB_MASK_3) << 1) |
		(find_0b110(opponent.b060, empty.b060) & BB_MASK_3);
	const uint64_t b060 = can_mate_1ply_1dir(player.b060, opponent.b060, empty.b060 & (~player_losing_060), opponent_losing_060);
	if (b060) {
		unsigned long pos = 0;
		bitscan_forward64(&pos, rotate_bitboard_300degree_clockwise(b060));
		bestmove = int(pos);
		return Immediate::WIN;
	}

	const uint64_t opponent_losing_120 =
		((find_0b011(opponent.b120, empty.b120) & BB_MASK_3) << 2) |
		((find_0b101(opponent.b120, empty.b120) & BB_MASK_3) << 1) |
		(find_0b110(opponent.b120, empty.b120) & BB_MASK_3);
	const uint64_t b120 = can_mate_1ply_1dir(player.b120, opponent.b120, empty.b120 & (~player_losing_120), opponent_losing_120);
	if (b120) {
		unsigned long pos = 0;
		bitscan_forward64(&pos, rotate_bitboard_240degree_clockwise(b120));
		bestmove = int(pos);
		return Immediate::WIN;
	}

	return Immediate::UNKNOWN;
}

//三手詰めかどうかの判定を行う。この局面はplayerの手番だとする。
//
//このターンに石を置いて即座に勝利できる場合、この関数は正確な値を返さない。exist_winning_pos_bitboard関数で先にチェックしておくこと。
//
//返り値がImmediate::WINである⇒playerの勝ちである。
//返り値がImmediate::LOSEである⇒playerの負けである。
int search_mate_3ply(Bitboard &player, Bitboard &opponent, Bitboard &empty, const uint64_t opponent_winning, uint64_t &player_losing_bitboard_000) noexcept {

	if (popcount64(empty.b000) >= 61)return Immediate::UNKNOWN;

	//どの空白マスに打っても3並びになるならば負け。あるいは、現在王手がかかっておらず、かつ一手詰めの王手が存在するなら勝ち。
	const int is_mate_1ply = search_mate_1ply(player, opponent, empty, opponent_winning, player_losing_bitboard_000);
	if (is_mate_1ply != Immediate::UNKNOWN)return is_mate_1ply;

	//
	//攻め方
	//

	const uint64_t player_check = get_check_pos_bitboard(player, empty);
	uint64_t bb_move = player_check & (~player_losing_bitboard_000);

	if (opponent_winning != 0) {//受け方から王手をかけられている場合
		bb_move &= opponent_winning;//逆王手だけを探索する。
	}

	for (unsigned long pos = 0; bitscan_forward64(&pos, bb_move); bb_move &= bb_move - 1) {

		assert(pos < 61);

		player.do_move(uint8_t(pos));
		empty.undo_move(uint8_t(pos));
		{
			//
			//受け方
			//

			//攻め方が王手をかけそこねたとき、ここにひっかかる。
			assert(exist_winning_pos_bitboard(opponent, empty) == 0);

			const uint64_t __opponent_winning = get_winning_pos_bitboard(player, empty);
			assert(__opponent_winning > 0);
			const uint64_t __player_losing = get_losing_pos_bitboard(opponent, empty);
			if (((__opponent_winning & __player_losing) == __opponent_winning)
				|| ((__opponent_winning & (__opponent_winning - 1)) != 0)) {
				player.undo_move(uint8_t(pos));
				empty.do_move(uint8_t(pos));
				return Immediate::WIN;
			}

			unsigned long __pos = 0;
			bitscan_forward64(&__pos, __opponent_winning);

			opponent.do_move(uint8_t(__pos));
			empty.undo_move(uint8_t(__pos));
			{
				//
				//攻め方
				//

				//受け方が王手を見逃したときのみここにひっかかる。二重王手はチェックしてある。
				assert(exist_winning_pos_bitboard(player, empty) == 0);

				//どの空白マスに打っても3並びになるならば負け。あるいは、現在王手がかかっておらず、かつ一手詰めの王手が存在するなら勝ち。
				const uint64_t ____opponent_winning = get_winning_pos_bitboard(opponent, empty);
				uint64_t ____player_losing_bitboard_000 = 0;
				const int ____is_mate_1ply = search_mate_1ply(player, opponent, empty, ____opponent_winning, ____player_losing_bitboard_000);
				if (____is_mate_1ply == Immediate::WIN) {
					opponent.undo_move(uint8_t(__pos));
					empty.do_move(uint8_t(__pos));
					player.undo_move(uint8_t(pos));
					empty.do_move(uint8_t(pos));
					return Immediate::WIN;
				}

			}
			opponent.undo_move(uint8_t(__pos));
			empty.do_move(uint8_t(__pos));

		}
		player.undo_move(uint8_t(pos));
		empty.do_move(uint8_t(pos));

	}

	return Immediate::UNKNOWN;
}

//詰め探索を行う。すなわち、攻め方（＝height=0のときの指し手）は王手だけを指すとして、受け方を詰ませられるかを探索する。
int search_mate(Bitboard &player, Bitboard &opponent, Bitboard &empty, int depth, int height = 0) noexcept {

	assert(height >= 0);

	if (empty.b000 == 0)return Immediate::DRAW;

	//rootで入力された局面が既に詰んでいる場合のみここにひっかかる。
	assert(is_immediate_end(player) == Immediate::UNKNOWN);

	//このターンに勝利できるならば返す。
	if (exist_winning_pos_bitboard(player, empty))return Immediate::WIN - height - 1;


	const uint64_t opponent_winning = get_winning_pos_bitboard(opponent, empty);

	const uint64_t player_losing_000 =
		((find_0b011(player.b000, empty.b000) & BB_MASK_3) << 2) |
		((find_0b101(player.b000, empty.b000) & BB_MASK_3) << 1) |
		(find_0b110(player.b000, empty.b000) & BB_MASK_3);
	const uint64_t player_losing_060 =
		((find_0b011(player.b060, empty.b060) & BB_MASK_3) << 2) |
		((find_0b101(player.b060, empty.b060) & BB_MASK_3) << 1) |
		(find_0b110(player.b060, empty.b060) & BB_MASK_3);
	const uint64_t player_losing_120 =
		((find_0b011(player.b120, empty.b120) & BB_MASK_3) << 2) |
		((find_0b101(player.b120, empty.b120) & BB_MASK_3) << 1) |
		(find_0b110(player.b120, empty.b120) & BB_MASK_3);
	const uint64_t player_losing = player_losing_000 | rotate_bitboard_300degree_clockwise(player_losing_060) | rotate_bitboard_240degree_clockwise(player_losing_120);

	if (depth <= 0)return Immediate::UNKNOWN;

	int bestscore = Immediate::LOSE;
	uint64_t bb_move = 0;

	if (height % 2 == 0) {//攻め方の指し手生成

		bestscore = Immediate::UNKNOWN;

		const uint64_t player_check = get_check_pos_bitboard(player, empty);
		uint64_t bb_move = player_check & (~player_losing);

		if (opponent_winning != 0) {//受け方から王手をかけられている場合
			bb_move &= opponent_winning;//逆王手だけを探索する。
		}
	}
	else {//受け方の合法手生成
		assert(opponent_winning != 0);//攻め方から必ず王手をかけられている。

		bb_move = opponent_winning;
	}

	for (unsigned long pos = 0; bitscan_forward64(&pos, bb_move); bb_move &= bb_move - 1) {

		assert(pos < 61);

		player.do_move(uint8_t(pos));
		empty.undo_move(uint8_t(pos));
		const int result = -search_mate(opponent, player, empty, depth - 1, height + 1);
		player.undo_move(uint8_t(pos));
		empty.do_move(uint8_t(pos));

		if (bestscore < result)bestscore = result;
	}

	return bestscore;
}

int search_endgame(Bitboard &player, Bitboard &opponent, Bitboard &empty, int alpha, const int beta) noexcept {

	if (empty.b000 == 0)return Immediate::DRAW;

	//rootで入力された局面が既に詰んでいる場合のみここにひっかかる。
	assert(is_immediate_end(player) == Immediate::UNKNOWN);

	//このターンに勝利できるならば返す。
	if (exist_winning_pos_bitboard(player, empty))return Immediate::WIN;

	const uint64_t opponent_winning = get_winning_pos_bitboard(opponent, empty);
	uint64_t player_losing = 0;
	const int is_mate_1ply = search_mate_1ply(player, opponent, empty, opponent_winning, player_losing);
	if (is_mate_1ply != Immediate::UNKNOWN)return is_mate_1ply;

	if (opponent_winning) {

		unsigned long pos = 0;
		bitscan_forward64(&pos, opponent_winning);

		player.do_move(uint8_t(pos));
		empty.undo_move(uint8_t(pos));
		const int result = -search_endgame(opponent, player, empty, -beta, -alpha);
		player.undo_move(uint8_t(pos));
		empty.do_move(uint8_t(pos));

		return result;
	}

	//空白マスのうち、打っても即負けにならない場所のビットボードを作る。
	uint64_t bb_move = empty.b000 & (~player_losing);

	//どの空白マスに打っても即負ける事態は既に検知済みなので大丈夫のはず。
	assert(bb_move > 0);

	int bestscore = Immediate::LOSE;

	for (unsigned long pos = 0; bitscan_forward64(&pos, bb_move); bb_move &= bb_move - 1) {

		assert(pos < 61);

		player.do_move(uint8_t(pos));
		empty.undo_move(uint8_t(pos));
		const int result = -search_endgame(opponent, player, empty, -beta, -alpha);
		player.undo_move(uint8_t(pos));
		empty.do_move(uint8_t(pos));

		if (bestscore < result)bestscore = result;
		if (alpha < result)alpha = result;
		if (beta <= result)return result;
	}

	return bestscore;
}

void search_endgame_root(Bitboard &player, Bitboard &opponent, Bitboard &empty, int &final_score, int& bestmove) {

	if (empty.b000 == 0) {
		final_score = Immediate::DRAW;
		bestmove = -1;
		return;
	}

	//このターンに勝利できるならば返す。
	if (exist_winning_pos_bitboard(player, empty)) {
		unsigned long pos = 0;
		bitscan_forward64(&pos, get_winning_pos_bitboard(player, empty));
		bestmove = int(pos);
		final_score = Immediate::WIN;
		return;
	}

	const uint64_t opponent_winning = get_winning_pos_bitboard(opponent, empty);

	int bestmove_1ply = 0;
	uint64_t player_losing = 0;
	const int is_mate_1ply = search_mate_1ply(player, opponent, empty, opponent_winning, player_losing, bestmove_1ply);
	if (is_mate_1ply != Immediate::UNKNOWN) {
		bestmove = bestmove_1ply;
		final_score = is_mate_1ply;
		return;
	}

	if (opponent_winning) {

		unsigned long pos = 0;
		bitscan_forward64(&pos, opponent_winning);

		player.do_move(uint8_t(pos));
		empty.undo_move(uint8_t(pos));
		const int result = -search_endgame(opponent, player, empty, Immediate::LOSE, Immediate::WIN);
		player.undo_move(uint8_t(pos));
		empty.do_move(uint8_t(pos));

		bestmove = int(pos);
		final_score = result;
		return;
	}

	//空白マスのうち、打っても即負けにならない場所のビットボードを作る。
	uint64_t bb_move = empty.b000 & (~player_losing);

	//どの空白マスに打っても即負ける事態は既に検知済みなので大丈夫のはず。
	assert(bb_move > 0);

	final_score = Immediate::LOSE;
	int alpha = Immediate::LOSE;

	for (unsigned long pos = 0; bitscan_forward64(&pos, bb_move); bb_move &= bb_move - 1) {

		assert(pos < 61);

		player.do_move(uint8_t(pos));
		empty.undo_move(uint8_t(pos));
		const int result = -search_endgame(opponent, player, empty, Immediate::LOSE, -alpha);
		player.undo_move(uint8_t(pos));
		empty.do_move(uint8_t(pos));

		if (final_score <= result) {
			final_score = result;
			bestmove = int(pos);
			if (alpha < result)alpha = result;
			if (result == Immediate::WIN) {
				final_score = Immediate::WIN;
				return;
			}
		}
	}

	return;
}

int search_midgame(Bitboard &player, Bitboard &opponent, Bitboard &empty, int depth) {

	const int immediate_player = is_immediate_end(player);
	const int immediate_opponent = is_immediate_end(opponent);

	if (immediate_player == Immediate::WIN) {
		return Immediate::WIN;
	}
	if (immediate_player == Immediate::LOSE) {
		return Immediate::LOSE;
	}
	if (immediate_opponent == Immediate::LOSE) {
		return Immediate::WIN;
	}
	if (immediate_opponent == Immediate::WIN) {
		return Immediate::LOSE;
	}

	if (empty.b000 == 0)return Immediate::DRAW;
	if (depth == 0)return Immediate::UNKNOWN;

	//このターンに勝利できるならば返す。
	if (exist_winning_pos_bitboard(player, empty))return Immediate::WIN;

	const uint64_t opponent_winning = get_winning_pos_bitboard(opponent, empty);
	uint64_t player_losing = 0;
	const int is_mate_1ply = search_mate_1ply(player, opponent, empty, opponent_winning, player_losing);
	if (is_mate_1ply != Immediate::UNKNOWN)return is_mate_1ply;

	if (opponent_winning) {

		unsigned long pos = 0;
		bitscan_forward64(&pos, opponent_winning);

		player.do_move(uint8_t(pos));
		empty.undo_move(uint8_t(pos));
		int result = -search_midgame(opponent, player, empty, depth - 1);
		player.undo_move(uint8_t(pos));
		empty.do_move(uint8_t(pos));

		if (result == -Immediate::UNKNOWN)return Immediate::UNKNOWN;
		return result;
	}

	//空白マスのうち、打っても即負けにならない場所のビットボードを作る。
	uint64_t bb_move = empty.b000 & (~player_losing);

	//どの空白マスに打っても即負ける事態は既に検知済みなので大丈夫のはず。
	assert(bb_move > 0);

	int bestscore = Immediate::LOSE;

	for (unsigned long pos = 0; bitscan_forward64(&pos, bb_move); bb_move &= bb_move - 1) {

		assert(pos < 61);

		player.do_move(uint8_t(pos));
		empty.undo_move(uint8_t(pos));
		int result = -search_midgame(opponent, player, empty, depth - 1);
		player.undo_move(uint8_t(pos));
		empty.do_move(uint8_t(pos));

		if (result == -Immediate::UNKNOWN)result = Immediate::UNKNOWN;

		if (bestscore < result)bestscore = result;
	}

	return bestscore;
}


//ミニマックス探索を行う。is_immediate_endが正しく動くことだけを仮定する。他の関数のテストのために使う。
int search_naive(Bitboard &player, Bitboard &opponent, Bitboard &empty, int depth, int height = 0) {

	const int immediate_player = is_immediate_end(player);
	const int immediate_opponent = is_immediate_end(opponent);

	if (immediate_player == Immediate::WIN) {
		return Immediate::WIN - height;
	}
	if (immediate_player == Immediate::LOSE) {
		return Immediate::LOSE + height;
	}
	if (immediate_opponent == Immediate::LOSE) {
		return Immediate::WIN - height;
	}
	if (immediate_opponent == Immediate::WIN) {
		return Immediate::LOSE + height;
	}

	if (empty.b000 == 0)return Immediate::DRAW;
	if (depth == 0)return Immediate::UNKNOWN;

	const uint64_t opponent_winning = get_winning_pos_bitboard(opponent, empty);

	uint64_t bb_move = empty.b000;

	int bestscore = Immediate::LOSE;

	if (opponent_winning) {

		//王手をかけられているときに、それを無視して石を置いて即勝利できるかどうかを調べる。
		//二重王手に対して逆王手で応じた場合に
		//
		//https://nestorgames.com/docs/YavalathCo2.pdf
		//メモ：3人向けルールでは、「自分の直後に指すプレイヤーが王手をかけているときは可能な限りそれを防がなければならない」というルールがある。
		//      ゆえに、その人の王手を無視して自分が即勝利することは許されない。
		//      以上は3人向けルール限定の事態であり、2人向けルールには関係ない（即勝利してよい）。
		//メモ：この関数の返り値の絶対値は、勝敗が確定するならば(1000-手数)である。つまり、勝てる側は最短で勝とうとして、負ける側は可能な限り粘ろうとする。
		//      王手を防ぐと3目並んで負けになる場合、2人向けyavalathのルールに厳密に従うなら、王手を無視して1手引き延ばすのが最善である。（そのような手が存在する場合に限る）
		//      しかし今回の実装では、王手を無視する手は（それで即勝利できる場合を限り）選"べ"ないと仮定して手数を計算する。

		uint64_t bb_move_ignoring_check = bb_move & (~opponent_winning);
		for (unsigned long pos = 0; bitscan_forward64(&pos, bb_move_ignoring_check); bb_move_ignoring_check &= bb_move_ignoring_check - 1) {

			assert(pos < 61);

			player.do_move(uint8_t(pos));
			empty.undo_move(uint8_t(pos));
			if (is_immediate_end(player) == Immediate::WIN) {
				const int result = Immediate::WIN - (height + 1);
				if (bestscore < result)bestscore = result;
			}
			player.undo_move(uint8_t(pos));
			empty.do_move(uint8_t(pos));
		}

		bb_move &= opponent_winning;
	}

	for (unsigned long pos = 0; bitscan_forward64(&pos, bb_move); bb_move &= bb_move - 1) {

		assert(pos < 61);

		player.do_move(uint8_t(pos));
		empty.undo_move(uint8_t(pos));
		const int result = -search_naive(opponent, player, empty, depth - 1, height + 1);
		player.undo_move(uint8_t(pos));
		empty.do_move(uint8_t(pos));

		if (bestscore < result)bestscore = result;
	}

	return bestscore;
}

int random_play_naive(Bitboard &player, Bitboard &opponent, Bitboard &empty, Xoroshiro128plus &rnd, std::array<int, 62> &kifu) {

	const int immediate_player = is_immediate_end(player);
	const int immediate_opponent = is_immediate_end(player);

	if (immediate_player == Immediate::WIN || is_immediate_end(opponent) == Immediate::LOSE)return Immediate::WIN;
	if (immediate_player == Immediate::LOSE || is_immediate_end(opponent) == Immediate::WIN)return Immediate::LOSE;

	if (empty.b000 == 0)return Immediate::DRAW;

	const uint64_t opponent_winning = get_winning_pos_bitboard(opponent, empty);

	if (opponent_winning) {

		unsigned long pos = 0;
		bitscan_forward64(&pos, opponent_winning);

		kifu[61 - popcount64(empty.b000)] = pos;
		player.do_move(uint8_t(pos));
		empty.undo_move(uint8_t(pos));
		return -random_play_naive(opponent, player, empty, rnd, kifu);
	}

	uint64_t bb_move = empty.b000;
	if (opponent_winning)bb_move &= opponent_winning;

	const uint64_t pop = popcount64(bb_move);

	const uint64_t order = rnd.rand() % pop;//order←0以上pop未満のランダムな値

	const uint64_t selected_bb_move = pdep64(1ULL << order, bb_move);
	unsigned long pos = 0;
	bitscan_forward64(&pos, selected_bb_move);

	kifu[61 - popcount64(empty.b000)] = pos;
	player.do_move(uint8_t(pos));
	empty.undo_move(uint8_t(pos));
	return -random_play_naive(opponent, player, empty, rnd, kifu);
}

int random_play_naive_root(uint64_t seed, std::array<int, 62> &kifu) {

	Xoroshiro128plus rnd(seed);

	Bitboard player, opponent, empty(BB_ALL);

	return random_play_naive(player, opponent, empty, rnd, kifu);
}

//引数局面からお互いランダムに指して、勝敗を返す。ただし完全ランダムではなく、即負けになる手や、直後に3手詰めされるような手は指さない。
std::pair<int, bool> random_play_pair(Bitboard &player, Bitboard &opponent, Bitboard &empty, Xoroshiro128plus &rnd) noexcept {

	//if (empty.b000 == 0)return Immediate::DRAW;
	if (popcount64(empty.b000) <= 6) {
		return std::make_pair(search_endgame(player, opponent, empty, Immediate::LOSE, Immediate::WIN), false);
	}


	//rootで入力された局面が既に詰んでいる場合のみここにひっかかる。
	assert(is_immediate_end(player) == Immediate::UNKNOWN);
	assert(is_immediate_end(opponent) == Immediate::UNKNOWN);

	//このターンに勝利できるならば返す。
	if (exist_winning_pos_bitboard(player, empty))return std::make_pair(Immediate::WIN, false);


	const uint64_t opponent_winning = get_winning_pos_bitboard(opponent, empty);
	uint64_t player_losing = 0;
	const int is_mate = search_mate_3ply(player, opponent, empty, opponent_winning, player_losing);
	if (is_mate == Immediate::LOSE)return std::make_pair(Immediate::LOSE, false);
	else if (is_mate == Immediate::WIN)return std::make_pair(Immediate::WIN, true);


	if (opponent_winning) {

		unsigned long pos = 0;
		bitscan_forward64(&pos, opponent_winning);

		player.do_move(uint8_t(pos));
		empty.undo_move(uint8_t(pos));

		const auto result = random_play_pair(opponent, player, empty, rnd);
		return std::make_pair(-result.first, false);
	}

	//空白マスのうち、打っても即負けにならない場所のビットボードを作る。
	uint64_t bb_move = empty.b000 & (~player_losing);

	//どの空白マスに打っても即負ける事態は既に検知済みなので大丈夫のはず。
	assert(bb_move > 0);


	//ランダムに指す。ただし、その直後に相手方の3手詰めで負けてしまう場合に限り、別の手をランダムに選び直す。
	for (uint64_t pop = popcount64(bb_move); pop; --pop) {

		const uint64_t order = rnd.rand() % pop;

		const uint64_t selected_bb_move = pdep64(1ULL << order, bb_move);
		assert(popcount64(selected_bb_move) == 1);
		unsigned long pos = 0;
		bitscan_forward64(&pos, selected_bb_move);

		player.do_move(uint8_t(pos));
		empty.undo_move(uint8_t(pos));
		const auto result = random_play_pair(opponent, player, empty, rnd);
		if (result.second == false)return std::make_pair(-result.first, false);
		player.undo_move(uint8_t(pos));
		empty.do_move(uint8_t(pos));

		assert(bb_move & selected_bb_move);
		bb_move ^= selected_bb_move;

	}

	//どの手を指しても相手方に3手詰めさせられてしまうので負け。
	return std::make_pair(Immediate::LOSE, false);
}



int random_play_root(uint64_t seed) {

	Xoroshiro128plus rnd(seed);

	Bitboard player, opponent, empty(BB_ALL);

	const auto result = random_play_pair(player, opponent, empty, rnd);

	return result.first;
}

//引数局面において、考慮する価値のある指し手のビットボードを返り値のfirstとする。引数局面の勝敗が明らかなら、firstはゼロとなる。
//返り値のsecondは、引数局面の勝敗が明らかならその勝敗で、さもなくばImmediate::UNKNOWN。
std::pair<uint64_t, int>get_all_reasonable_moves(Bitboard &player, Bitboard &opponent, Bitboard &empty) {

	if (popcount64(empty.b000) <= 6) {
		return std::make_pair(0, search_endgame(player, opponent, empty, Immediate::LOSE, Immediate::WIN));
	}

	//rootで入力された局面が既に詰んでいる場合のみここにひっかかる。
	assert(is_immediate_end(player) == Immediate::UNKNOWN);
	assert(is_immediate_end(opponent) == Immediate::UNKNOWN);

	//このターンに勝利できるならば返す。
	if (exist_winning_pos_bitboard(player, empty))return std::make_pair(0, Immediate::WIN);

	const uint64_t opponent_winning = get_winning_pos_bitboard(opponent, empty);
	uint64_t player_losing = 0;
	const int is_mate = search_mate_3ply(player, opponent, empty, opponent_winning, player_losing);
	if (is_mate != Immediate::UNKNOWN)return std::make_pair(0, is_mate);

	uint64_t bb_move = empty.b000 & (~player_losing);
	if (opponent_winning)bb_move = opponent_winning;

	//どの空白マスに打っても即負ける事態は既に検知済みなので大丈夫のはず。
	assert(bb_move > 0);

	uint64_t answer = 0;

	for (unsigned long pos = 0; bitscan_forward64(&pos, bb_move); bb_move &= bb_move - 1) {

		player.do_move(uint8_t(pos));
		empty.undo_move(uint8_t(pos));
		const auto result = -search_midgame(opponent, player, empty, 1);
		player.undo_move(uint8_t(pos));
		empty.do_move(uint8_t(pos));

		if (result == Immediate::WIN)return std::make_pair(0, Immediate::WIN);
		if (result == Immediate::LOSE)continue;

		answer |= 1ULL << pos;
	}

	if (answer == 0) {
		return std::make_pair(0, Immediate::LOSE);
	}

	return std::make_pair(answer, Immediate::UNKNOWN);
}

struct MCTS_Node {
	uint64_t player_b000, opponent_b000, reasonable_moves_b000, unvisited_moves_b000;
	int64_t visit_count, win_count, lose_count;
	uint64_t next_start_index, parent_index;
	int confirmed_result;

	MCTS_Node(const uint64_t player, const uint64_t opponent) {
		player_b000 = player;
		opponent_b000 = opponent;
		reasonable_moves_b000 = unvisited_moves_b000 = 0;
		visit_count = win_count = lose_count = 0;
		next_start_index = parent_index = -1;
		confirmed_result = Immediate::UNKNOWN;
	}

	MCTS_Node() :MCTS_Node(0, 0) {}
};

constexpr double SQUARE_OF_UCT_CONST_C = 0.5;

bool mcts_stop_flag = false;

template<bool output_mind>class MCTS_Solver {

	const uint64_t Nodes_limit = 10000000;

	std::vector<MCTS_Node>Nodes;

	Xoroshiro128plus rnd;

	static int answer_if_obvious_position(Bitboard player, Bitboard opponent, Bitboard empty) {

		//引数局面について、最善手が自明な局面だった場合にはその手(0以上61未満)を返す。さもなくば-1を返す。ただし、指しようがない局面に対しては例外を投げる。
		//自明な局面とは：
		//(1)空白が6箇所以下の場合、厳密解を求めて返す。
		//(2)1手詰めが検出された場合、厳密解を求めて返す。
		//(3)王手をかけられていて合法手が1手しか無い場合、その手を返す。
		//(4)各合法手の先の局面を1手読みした結果、
		//(4-1)勝てる手が見つかった場合、最初に見つかった手を返す。
		//(4-2)どの手を指しても負けるとわかった場合、即負けにならない手のうちどれかを返す。
		//(4-3)負け確定だと言い切れない手が1つだけある場合、その手を返す。
		//以上のこととする。


		const auto num_empty = popcount64(empty.b000);

		if (num_empty == 0 ||
			is_immediate_end(player) != Immediate::UNKNOWN ||
			is_immediate_end(opponent) != Immediate::UNKNOWN) {
			throw std::runtime_error("think_at_first: root node has already finished");
		}

		if (popcount64(empty.b000) <= 6) {
			int finalscore = 0, bestmove = 0;
			search_endgame_root(player, opponent, empty, finalscore, bestmove);

			if (output_mind) {
				if (finalscore == Immediate::WIN)std::cout << "obvious_move: endgame is solved: computer is confident of victory." << std::endl;
				else if (finalscore == Immediate::LOSE)std::cout << "obvious_move: endgame is solved: computer is predicting its defeat." << std::endl;
				else if (finalscore == Immediate::DRAW)std::cout << "obvious_move: endgame is solved: computer is predicting a draw." << std::endl;
				else assert(false);
			}

			return bestmove;
		}

		//このターンに勝利できるならば返す。
		if (exist_winning_pos_bitboard(player, empty)) {
			if (output_mind)std::cout << "obvious_move: winning move: computer is confident of victory." << std::endl;
			unsigned long pos = 0;
			bitscan_forward64(&pos, get_winning_pos_bitboard(player, empty));
			return int(pos);
		}

		const uint64_t opponent_winning = get_winning_pos_bitboard(opponent, empty);

		int bestmove_1ply = 0;
		uint64_t player_losing = 0;
		const int is_mate_1ply = search_mate_1ply(player, opponent, empty, opponent_winning, player_losing, bestmove_1ply);
		if (is_mate_1ply != Immediate::UNKNOWN) {
			if (output_mind) {
				if (is_mate_1ply == Immediate::WIN)std::cout << "obvious_move: search_mate_1ply: computer is confident of victory." << std::endl;
				else if (is_mate_1ply == Immediate::LOSE)std::cout << "obvious_move: search_mate_1ply: computer is predicting its defeat." << std::endl;
				else assert(false);
			}

			return bestmove_1ply;
		}

		if (opponent_winning) {
			assert(popcount64(opponent_winning) == 1);
			unsigned long pos = 0;
			bitscan_forward64(&pos, opponent_winning);
			if (output_mind)std::cout << "obvious_move: forced reply." << std::endl;
			return int(pos);
		}

		uint64_t bb_move = empty.b000 & (~player_losing);

		//どの空白マスに打っても即負ける事態は既に検知済みなので大丈夫のはず。
		assert(bb_move > 0);

		uint64_t answer = 0;

		for (unsigned long pos = 0; bitscan_forward64(&pos, bb_move); bb_move &= bb_move - 1) {

			player.do_move(uint8_t(pos));
			empty.undo_move(uint8_t(pos));
			const auto result = -search_midgame(opponent, player, empty, 1);
			player.undo_move(uint8_t(pos));
			empty.do_move(uint8_t(pos));

			if (result == Immediate::WIN) {
				if (output_mind)std::cout << "obvious_move: search_midgame(depth=1): computer is confident of victory." << std::endl;
				return int(pos);
			}
			if (result == Immediate::LOSE)continue;

			answer |= 1ULL << pos;
		}

		if (answer == 0) {
			unsigned long pos = 0;
			bitscan_forward64(&pos, empty.b000 & (~player_losing));
			if (output_mind)std::cout << "obvious_move: search_midgame(depth=1): computer is predicting its defeat." << std::endl;
			return int(pos);
		}
		if (popcount64(answer) == 1) {
			unsigned long pos = 0;
			bitscan_forward64(&pos, answer);
			if (output_mind)std::cout << "obvious_move: search_midgame(depth=1): computer is predicting that all but one move(s) will result in defeat." << std::endl;
			return int(pos);
		}

		return -1;
	}

public:

	MCTS_Solver(uint64_t root_player_b000, uint64_t root_opponent_b000, uint64_t random_seed) :rnd(random_seed) {

		Nodes = std::vector<MCTS_Node>{ MCTS_Node(root_player_b000, root_opponent_b000) };

		Nodes.reserve(Nodes_limit + 100);

	}

private:

	uint64_t selection(const uint64_t node_index) {
		//Nodes[node_index]が葉ではないと仮定する。
		//子ノードの勝敗確定を見てNodes[node_index]の勝敗が確定する場合はnode_indexを返す。
		//さもなくば、勝敗が確定していない子ノードを以下の基準により選んでindexを返す。
		//(1)Nodes[node_index]の子のなかに葉があるなら、葉のどれかを選ぶ。
		//(2)Nodes[node_index]の子のなかに葉が無いなら、UCT値が最も大きい子を選ぶ。

		assert(Nodes[node_index].visit_count > 0);

		if (Nodes[node_index].unvisited_moves_b000 != 0) {
			const uint64_t pop = popcount64(Nodes[node_index].unvisited_moves_b000);
			const uint64_t order_among_unvisited = rnd.rand() % pop;
			const uint64_t selected_bb_move = pdep64(1ULL << order_among_unvisited, Nodes[node_index].unvisited_moves_b000);

			unsigned long order_among_all = 0;
			bitscan_forward64(&order_among_all, pext64(selected_bb_move, Nodes[node_index].reasonable_moves_b000));

			Nodes[node_index].unvisited_moves_b000 ^= selected_bb_move;
			return Nodes[node_index].next_start_index + order_among_all;
		}

		double UCT_max = std::numeric_limits<double>::lowest();
		int index_max = -1;
		int confirmed_best_result = Immediate::LOSE;

		const uint64_t pop = popcount64(Nodes[node_index].reasonable_moves_b000);
		const double C_C_log2_N = SQUARE_OF_UCT_CONST_C * std::log2(Nodes[node_index].visit_count);
		for (int i = 0; i < pop; ++i) {
			const uint64_t i_th_move_index = Nodes[node_index].next_start_index + i;

			if (Nodes[i_th_move_index].confirmed_result != Immediate::UNKNOWN) {
				if (Nodes[i_th_move_index].confirmed_result == Immediate::WIN) {
					continue;
				}
				else if (Nodes[i_th_move_index].confirmed_result == Immediate::LOSE) {
					Nodes[node_index].confirmed_result = Immediate::WIN;
					return node_index;
				}
				else if (Nodes[i_th_move_index].confirmed_result == Immediate::DRAW) {
					confirmed_best_result = Immediate::DRAW;
					continue;
				}
				else {
					assert(false);
				}
			}

			//勝利は+1.0, 引き分けは+0.5として換算する。式変形すると以下のように書ける。
			const double w_i = (Nodes[i_th_move_index].visit_count + Nodes[i_th_move_index].win_count - Nodes[i_th_move_index].lose_count) / 2.0;

			const double n_i = Nodes[i_th_move_index].visit_count;

			const double UCT_value = (w_i / n_i) + std::sqrt(C_C_log2_N / n_i);

			if (UCT_max < UCT_value) {
				UCT_max = UCT_value;
				index_max = i;
			}
		}

		//どの手を指してもその直後に相手の勝ちor引き分けが確定するのなら、この局面で自分の負けor引き分けが確定している。
		if (index_max == -1) {
			Nodes[node_index].confirmed_result = confirmed_best_result;
			return node_index;
		}

		return Nodes[node_index].next_start_index + index_max;
	}

	uint64_t expansion(const uint64_t node_index) {
		//Nodes[node_index]が葉だと仮定する。
		//Nodes[node_index]の勝敗が確定しているか調べ、いるならその勝敗を記録する。
		//さもなくば、合法手を全列挙して子ノードを生成する。これによりNodes[node_index]は葉でなくなる。
		//合法手が一つしかない場合（e.g.,王手逃れ）、更にexpansionする。
		//返り値は最終的に到達した『複数の子を持ち、かつそれら全てが葉であるような、自らは葉でないノード』(もしくは勝敗が確定したノード)のindexとする。

		assert(Nodes[node_index].reasonable_moves_b000 == 0);
		assert(Nodes[node_index].unvisited_moves_b000 == 0);
		assert(Nodes[node_index].visit_count == 0);
		assert(Nodes[node_index].win_count == 0);
		assert(Nodes[node_index].lose_count == 0);
		assert(Nodes[node_index].next_start_index == -1);
		assert(Nodes[node_index].parent_index != -1 || node_index == 0);

		if (Nodes[node_index].confirmed_result != Immediate::UNKNOWN)return node_index;

		Bitboard player(Nodes[node_index].player_b000);
		Bitboard opponent(Nodes[node_index].opponent_b000);
		Bitboard empty((~(Nodes[node_index].player_b000 | Nodes[node_index].opponent_b000)) & BB_ALL);

		//返り値のfirstは合法手のビットボード。secondは勝敗が確定しているかどうかの情報。
		const std::pair<uint64_t, int> moves = get_all_reasonable_moves(player, opponent, empty);

		//勝敗が確定しているならそれを記録して終了。
		if (moves.second != Immediate::UNKNOWN) {
			Nodes[node_index].confirmed_result = moves.second;
			return node_index;
		}

		assert(moves.first != 0);

		Nodes[node_index].reasonable_moves_b000 = moves.first;
		Nodes[node_index].unvisited_moves_b000 = moves.first;
		Nodes[node_index].next_start_index = Nodes.size();

		//各合法手を指した先の局面を子ノードとして登録する。
		uint64_t m = moves.first;
		for (unsigned long pos = 0; bitscan_forward64(&pos, m); m &= m - 1) {
			player.do_move(uint8_t(pos));
			Nodes.push_back(MCTS_Node{ opponent.b000, player.b000 });
			Nodes.back().parent_index = node_index;
			player.undo_move(uint8_t(pos));
		}

		if (popcount64(moves.first) == 1)return expansion(Nodes[node_index].next_start_index);

		return node_index;
	}

	int simulation(const uint64_t node_index) {

		assert(Nodes[node_index].reasonable_moves_b000 != 0);
		assert(Nodes[node_index].unvisited_moves_b000 != 0);
		assert(Nodes[node_index].visit_count == 0);
		assert(Nodes[node_index].win_count == 0);
		assert(Nodes[node_index].lose_count == 0);
		assert(Nodes[node_index].next_start_index != -1);

		Bitboard player(Nodes[node_index].player_b000);
		Bitboard opponent(Nodes[node_index].opponent_b000);
		Bitboard empty((~(Nodes[node_index].player_b000 | Nodes[node_index].opponent_b000)) & BB_ALL);

		return random_play_pair(player, opponent, empty, rnd).first;
	}

	void backpropagation_simulation(int result, uint64_t node_index) {

		while (node_index != -1) {

			++Nodes[node_index].visit_count;

			if (result > 0) {
				++Nodes[node_index].win_count;
			}
			else if (result < 0) {
				++Nodes[node_index].lose_count;
			}

			node_index = Nodes[node_index].parent_index;
			result *= -1;
		}
	}

	void backpropagation_confirmed(uint64_t node_index) {

		static_assert(Immediate::LOSE == -Immediate::WIN, "static_assert: Immediate::LOSE == -Immediate::WIN");
		static_assert(Immediate::DRAW == 0, "static_assert: Immediate::DRAW == 0");
		static_assert(Immediate::DRAW != Immediate::WIN, "static_assert: Immediate enums are different each other");
		static_assert(Immediate::DRAW != Immediate::LOSE, "static_assert: Immediate enums are different each other");
		static_assert(Immediate::DRAW != Immediate::UNKNOWN, "static_assert: Immediate enums are different each other");
		static_assert(Immediate::WIN != Immediate::LOSE, "static_assert: Immediate enums are different each other");
		static_assert(Immediate::WIN != Immediate::UNKNOWN, "static_assert: Immediate enums are different each other");
		static_assert(Immediate::LOSE != Immediate::UNKNOWN, "static_assert: Immediate enums are different each other");

		int result = Nodes[node_index].confirmed_result;
		assert(result != Immediate::UNKNOWN);

		node_index = Nodes[node_index].parent_index;
		result *= -1;

		while (node_index != -1) {

			if (result == Immediate::WIN) {
				assert(
					Nodes[node_index].confirmed_result == Immediate::WIN ||
					Nodes[node_index].confirmed_result == Immediate::UNKNOWN
				);
				Nodes[node_index].confirmed_result = Immediate::WIN;
			}
			else if (popcount64(Nodes[node_index].reasonable_moves_b000) == 1) {
				Nodes[node_index].confirmed_result = result;
			}
			else {
				break;
				//複数の合法手があったが、今回の手以外の全ての手の勝敗が既に確定していて、
				//今回の手も確定したことによりこのノードの勝敗も確定するというケースがある。けれども、
				//その検証は重いのでselection関数で行う。これはある種の遅延評価といえる。
			}

			node_index = Nodes[node_index].parent_index;
			result *= -1;
		}
	}

	void MCTS_one_rollout() {

		uint64_t node_index = 0;//root

		while (Nodes[node_index].visit_count > 0 && Nodes[node_index].confirmed_result == Immediate::UNKNOWN) {
			node_index = selection(node_index);
		}

		if (Nodes[node_index].confirmed_result == Immediate::UNKNOWN) {
			expansion(node_index);
		}

		if (Nodes[node_index].confirmed_result == Immediate::UNKNOWN) {
			const int result = simulation(node_index);
			backpropagation_simulation(result, node_index);
		}
		else {
			backpropagation_confirmed(node_index);
		}
	}

	int get_the_best_move() {

		int best_move = -1;
		double best_score = std::numeric_limits<double>::lowest();

		uint64_t bb_move = Nodes[0].reasonable_moves_b000;
		for (unsigned long pos = 0, count = 0; bitscan_forward64(&pos, bb_move); bb_move &= bb_move - 1, ++count) {
			const uint64_t i_th_move_index = Nodes[0].next_start_index + count;
			if (Nodes[i_th_move_index].confirmed_result == Immediate::LOSE) {
				if(output_mind)std::cout << "computer is confident of victory." << std::endl;
				return pos;
			}

			//勝利は+1.0, 引き分けは+0.5として換算する。式変形すると以下のように書ける。
			const double w_i = (Nodes[i_th_move_index].visit_count + Nodes[i_th_move_index].win_count - Nodes[i_th_move_index].lose_count) / 2.0;

			const double n_i = Nodes[i_th_move_index].visit_count;

			double move_score = 0.0;
			if (Nodes[i_th_move_index].confirmed_result == Immediate::DRAW)move_score = 0.5;
			else if (Nodes[i_th_move_index].confirmed_result == Immediate::WIN)move_score = -1.0 / (n_i + 1.0);//負けが確定している手の場合、多く訪れた局面のほうが複雑なのでマシだと言えるはず。
			else if (n_i > 0)move_score = w_i / n_i;

			if (best_score < move_score) {
				best_score = move_score;
				best_move = pos;
			}
		}

		assert(best_move != -1);

		if (output_mind)std::cout << "score = " << best_score << std::endl;

		return best_move;
	}

public:

	int MCTS_think_root(const int64_t limit_millisecond) {

		//自明な局面かどうか調べる。そうなら結果を求めて返す。
		{
			Bitboard player(Nodes[0].player_b000);
			Bitboard opponent(Nodes[0].opponent_b000);
			Bitboard empty((~(Nodes[0].player_b000 | Nodes[0].opponent_b000)) & BB_ALL);
			const int obvious_move = answer_if_obvious_position(player, opponent, empty);
			if (obvious_move != -1) {
				return obvious_move;
			}
		}

		//rolloutをやって、停止指示が来てないか確認するというのを繰り返す。たまに現状を標準出力する。
		const auto start = std::chrono::system_clock::now();
		auto start_rap = std::chrono::system_clock::now();
		int64_t count = 1;
		for (;; ++count) {
			MCTS_one_rollout();
			if (mcts_stop_flag)break;
			if (Nodes[0].confirmed_result != Immediate::UNKNOWN)break;
			if (count % 2 == 0) {
				const auto now = std::chrono::system_clock::now();

				const auto elapsed_time_rap = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_rap).count();
				if (elapsed_time_rap >= 1000) {
					start_rap = std::chrono::system_clock::now();


					//TODO:1秒おきに現状を標準出力で知らせる。


				}

				//引数の制限時間に達したら終了。
				const auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
				if (elapsed_time >= limit_millisecond) {
					break;
				}

				if (Nodes.size() >= Nodes_limit) {
					break;
				}

			}
		}

		if (output_mind)std::cout << count << " simulations in " << (double(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start).count()) / 1000.0) << "second" << std::endl;

		//最終結果を求めて返す。
		return get_the_best_move();
	}

};

void CLI_play_game(int64_t time_limit_millisecond) {

	Bitboard player, opponent, empty(BB_ALL);

	const auto is_end = [&]() {
		const int immediate_player = is_immediate_end(player);
		const int immediate_opponent = is_immediate_end(player);

		if (immediate_player == Immediate::WIN || is_immediate_end(opponent) == Immediate::LOSE) {
			std::cout << "player win!" << std::endl;
			print_board(player.b000, opponent.b000);
			return true;
		}
		if (immediate_player == Immediate::LOSE || is_immediate_end(opponent) == Immediate::WIN) {
			std::cout << "computer win!" << std::endl;
			print_board(player.b000, opponent.b000);
			return true;
		}
		if (popcount64(empty.b000) == 0) {
			std::cout << "draw!" << std::endl;
			print_board(player.b000, opponent.b000);
			return true;
		}
		return false;
	};

	while (true) {

		if (is_end())break;

		std::cout << "your turn" << std::endl;
		print_board(player.b000, opponent.b000);

		while (true) {
			std::string s;
			std::cin >> s;

			int pos = 0;
			try {
				pos = str2pos_bitboard(s);
				if ((empty.b000 & (1ULL << pos)) == 0) {
					std::cout << "error: occupied position. please retry." << std::endl;
				}
				player.do_move(uint8_t(pos));
				empty.undo_move(uint8_t(pos));
				break;
			}
			catch (...) {
				std::cout << "error: invalid string. please retry." << std::endl;
			}
		}

		if (is_end())break;

		std::cout << "computer's turn" << std::endl;
		print_board(player.b000, opponent.b000);

		{
			MCTS_Solver mcts = MCTS_Solver<true>(opponent.b000, player.b000, 12345);
			const int pos = mcts.MCTS_think_root(time_limit_millisecond);

			assert((empty.b000 & (1ULL << pos)) != 0);

			opponent.do_move(uint8_t(pos));
			empty.undo_move(uint8_t(pos));
		}
	}

	std::cout << "this program will exit after 60 seconds..." << std::endl;

	try {
		std::this_thread::sleep_for(std::chrono::seconds(60));
	}
	catch (...) {
		return;
	}
}

int mcts_vs_random(const int64_t time_limit_millisecond, const uint64_t random_seed, const bool random_moves_at_first) {
	//対局する。
	//mctsが勝ったら1を返し、randomが勝ったら-1を返し、引き分けなら0を返す。

	Bitboard player, opponent, empty(BB_ALL);

	Xoroshiro128plus rnd(random_seed);

	const auto is_end = [&]() {
		const int immediate_player = is_immediate_end(player);
		const int immediate_opponent = is_immediate_end(player);

		if (immediate_player == Immediate::WIN || is_immediate_end(opponent) == Immediate::LOSE) {
			std::cout << "random win!" << std::endl;
			print_board(player.b000, opponent.b000);
			return std::make_pair(true, -1);
		}
		if (immediate_player == Immediate::LOSE || is_immediate_end(opponent) == Immediate::WIN) {
			std::cout << "mcts win!" << std::endl;
			print_board(player.b000, opponent.b000);
			return std::make_pair(true, 1);
		}
		if (popcount64(empty.b000) == 0) {
			std::cout << "draw!" << std::endl;
			print_board(player.b000, opponent.b000);
			return std::make_pair(true, 0);
		}
		return std::make_pair(false, 0);
	};

	for(int i = 0;;++i) {

		if (i != 0 || random_moves_at_first) {

			{
				const auto e = is_end();
				if (e.first)return e.second;
			}

			std::cout << "random's turn" << std::endl;
			print_board(player.b000, opponent.b000);

			{
				const auto x = get_all_reasonable_moves(player, opponent, empty);
				if (x.second == Immediate::WIN) {
					std::cout << "random win!" << std::endl;
					return -1;
				}
				if (x.second == Immediate::LOSE) {
					std::cout << "mcts win!" << std::endl;
					return 1;
				}
				if (x.second == Immediate::DRAW) {
					std::cout << "draw!" << std::endl;
					return 0;
				}

				const uint64_t order = rnd.rand() % popcount64(x.first);

				const uint64_t selected_bb_move = pdep64(1ULL << order, x.first);
				assert(popcount64(selected_bb_move) == 1);
				unsigned long pos = 0;
				bitscan_forward64(&pos, selected_bb_move);

				player.do_move(uint8_t(pos));
				empty.undo_move(uint8_t(pos));
			}
		}


		{
			const auto e = is_end();
			if (e.first)return e.second;
		}

		std::cout << "mcts's turn" << std::endl;
		print_board(player.b000, opponent.b000);

		{
			MCTS_Solver mcts = MCTS_Solver<true>(opponent.b000, player.b000, rnd.rand());
			const int pos = mcts.MCTS_think_root(time_limit_millisecond);

			assert((empty.b000 & (1ULL << pos)) != 0);

			opponent.do_move(uint8_t(pos));
			empty.undo_move(uint8_t(pos));
		}
	}

	assert(false);
	return 0;
}

void test_mcts_strength(int iteration, uint64_t mcts_timelimit_milliseconds) {
	std::cout << "start: " << iteration << " random vs mcts from the initial position" << std::endl;
	const auto start = std::chrono::system_clock::now();
	int x[3] = {}, y[3] = {};
	for (int i = 0; i < iteration; ++i) {
		x[mcts_vs_random(mcts_timelimit_milliseconds, i, true) + 1]++;
		y[mcts_vs_random(mcts_timelimit_milliseconds, i + iteration, false) + 1]++;
	}
	const auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start).count();
	std::cout << "elapsed time: " << elapsed_time << " ms" << std::endl;
	std::cout << "black random win(s): " << x[0] << std::endl;
	std::cout << "white mcts win(s): " << x[2] << std::endl;
	std::cout << "draw(s)     : " << x[1] << std::endl;
	std::cout << "black mcts win(s): " << y[2] << std::endl;
	std::cout << "white random win(s): " << y[0] << std::endl;
	std::cout << "draw(s)     : " << y[1] << std::endl;

}

void test_search_functions(uint64_t seed, const int test_num = 100) {

	std::cout << "start: test_search_functions" << std::endl;

	auto onetest = [](Bitboard &player, Bitboard &opponent, Bitboard &empty) {

		if (empty.b000 == 0)return true;
		if (is_immediate_end(player) != Immediate::UNKNOWN)return true;
		if (is_immediate_end(opponent) != Immediate::UNKNOWN)return true;

		if (exist_winning_pos_bitboard(player, empty)) {
			const int result = search_naive(player, opponent, empty, 2);
			if (result < Immediate::WIN - 10) {
				std::cout << "error: exist_winning_pos_bitboard is positive / search_naive is negative" << std::endl;
				std::cout << "search result = " << result << std::endl;
				print_board(player.b000, opponent.b000);
				return false;
			}
		}
		else {
			const uint64_t opponent_winning = get_winning_pos_bitboard(opponent, empty);
			{
				uint64_t player_losing = 0;
				const int is_mate_1ply = search_mate_1ply(player, opponent, empty, opponent_winning, player_losing);
				if (is_mate_1ply != Immediate::UNKNOWN) {
					const int search_naive_result = search_naive(player, opponent, empty, 2);
					const int result1 = std::clamp(search_naive_result, -1, 1) * 1000;
					if (is_mate_1ply != result1) {
						std::cout << "error: search_mate_1ply != search_naive" << std::endl;
						std::cout << "search_mate_1ply result = " << is_mate_1ply << std::endl;
						std::cout << "search_naive result = " << search_naive_result << std::endl;
						std::cout << "result1 = " << result1 << std::endl;
						print_board(player.b000, opponent.b000);

						search_mate_1ply(player, opponent, empty, opponent_winning, player_losing);
						search_naive(player, opponent, empty, 2);

						return false;
					}
				}
			}
			{
				uint64_t player_losing = 0;
				const int is_mate_3ply = search_mate_3ply(player, opponent, empty, opponent_winning, player_losing);
				if (is_mate_3ply != Immediate::UNKNOWN) {
					const int search_naive_result = search_naive(player, opponent, empty, 4);
					const int result1 = std::clamp(search_naive_result, -1, 1) * 1000;
					if (is_mate_3ply != result1) {
						std::cout << "error: search_mate_3ply != search_naive" << std::endl;
						std::cout << "search_mate_3ply result = " << is_mate_3ply << std::endl;
						std::cout << "search_naive result = " << search_naive_result << std::endl;
						std::cout << "result1 = " << result1 << std::endl;
						print_board(player.b000, opponent.b000);

						search_mate_3ply(player, opponent, empty, opponent_winning, player_losing);
						search_naive(player, opponent, empty, 4);

						return false;
					}
				}
			}

		}

		return true;
	};

	std::array<int, 62>kifu;
	std::mt19937_64 rnd(seed);
	for (int i = 0; i < test_num; ++i) {

		//std::cout << i << std::endl;

		for (int j = 0; j < 62; ++j)kifu[j] = -1;

		random_play_naive_root(rnd(), kifu);

		Bitboard player1, player2, empty(BB_ALL);

		for (int j = 0; kifu[j] != -1; ++j) {

			if (j % 2 == 0) {
				player1.do_move(kifu[j]);
			}
			else {
				player2.do_move(kifu[j]);
			}

			empty.undo_move(kifu[j]);

			if (j % 2 == 0) {
				if (!onetest(player2, player1, empty))return;
			}
			else {
				if (!onetest(player1, player2, empty))return;
			}
		}
	}



	std::cout << "clear: test_search_functions" << std::endl;

}

void test_find_functions(uint64_t seed, const int test_num = 1000000) {

	std::cout << "start: test_find_functions" << std::endl;

	auto onetest = [](uint64_t x, uint64_t y) {
		{
			const uint64_t ans1 = find_bitpattern_naive(x, 0, 0b1111, 4);
			const uint64_t ans2 = find_4_or_more_consecutive_bits(x);
			if (ans1 != ans2) {
				std::cout << "error: find_4_or_more_consecutive_bits" << std::endl;
				std::cout << "x    = " << std::bitset<64>(x) << std::endl;
				std::cout << "ans1 = " << std::bitset<64>(ans1) << std::endl;
				std::cout << "ans2 = " << std::bitset<64>(ans2) << std::endl;
				return false;
			}
		}
		{
			const uint64_t ans1 = find_bitpattern_naive(x, 0, 0b111, 3);
			const uint64_t ans2 = find_3_or_more_consecutive_bits(x);
			if (ans1 != ans2) {
				std::cout << "error: find_3_or_more_consecutive_bits" << std::endl;
				std::cout << "x    = " << std::bitset<64>(x) << std::endl;
				std::cout << "ans1 = " << std::bitset<64>(ans1) << std::endl;
				std::cout << "ans2 = " << std::bitset<64>(ans2) << std::endl;
				return false;
			}
		}
		{
			const uint64_t ans1 = find_0b1011(x, y);
			const uint64_t ans2 = find_bitpattern_naive(x, y, 0b1011, 4);
			if (ans1 != ans2) {
				std::cout << "error: find_0b1011_naive" << std::endl;
				std::cout << "x    = " << std::bitset<64>(x) << std::endl;
				std::cout << "y    = " << std::bitset<64>(y) << std::endl;
				std::cout << "ans1 = " << std::bitset<64>(ans1) << std::endl;
				std::cout << "ans2 = " << std::bitset<64>(ans2) << std::endl;
				return false;
			}
		}
		{
			const uint64_t ans1 = find_0b1101(x, y);
			const uint64_t ans2 = find_bitpattern_naive(x, y, 0b1101, 4);
			if (ans1 != ans2) {
				std::cout << "error: find_0b1101_naive" << std::endl;
				std::cout << "x    = " << std::bitset<64>(x) << std::endl;
				std::cout << "y    = " << std::bitset<64>(y) << std::endl;
				std::cout << "ans1 = " << std::bitset<64>(ans1) << std::endl;
				std::cout << "ans2 = " << std::bitset<64>(ans2) << std::endl;
				return false;
			}
		}
		{
			const uint64_t ans1 = find_0b101(x, y);
			const uint64_t ans2 = find_bitpattern_naive(x, y, 0b101, 3);
			if (ans1 != ans2) {
				std::cout << "error: find_0b101_naive" << std::endl;
				std::cout << "x    = " << std::bitset<64>(x) << std::endl;
				std::cout << "y    = " << std::bitset<64>(y) << std::endl;
				std::cout << "ans1 = " << std::bitset<64>(ans1) << std::endl;
				std::cout << "ans2 = " << std::bitset<64>(ans2) << std::endl;
				return false;
			}
		}
		{
			const uint64_t ans1 = find_0b110(x, y);
			const uint64_t ans2 = find_bitpattern_naive(x, y, 0b110, 3);
			if (ans1 != ans2) {
				std::cout << "error: find_0b110_naive" << std::endl;
				std::cout << "x    = " << std::bitset<64>(x) << std::endl;
				std::cout << "y    = " << std::bitset<64>(y) << std::endl;
				std::cout << "ans1 = " << std::bitset<64>(ans1) << std::endl;
				std::cout << "ans2 = " << std::bitset<64>(ans2) << std::endl;
				return false;
			}
		}
		{
			const uint64_t ans1 = find_0b011(x, y);
			const uint64_t ans2 = find_bitpattern_naive(x, y, 0b011, 3);
			if (ans1 != ans2) {
				std::cout << "error: find_0b011_naive" << std::endl;
				std::cout << "x    = " << std::bitset<64>(x) << std::endl;
				std::cout << "y    = " << std::bitset<64>(y) << std::endl;
				std::cout << "ans1 = " << std::bitset<64>(ans1) << std::endl;
				std::cout << "ans2 = " << std::bitset<64>(ans2) << std::endl;
				return false;
			}
		}
		return true;
	};

	std::mt19937_64 rnd(seed);
	for (int i = 0; i < test_num; ++i) {

		{
			const uint64_t x = rnd();
			const uint64_t y = rnd() & (~x);
			assert((x & y) == 0);
			if (!onetest(x, y))return;
			if (!onetest(x, 0))return;
			if (!onetest(x, ~x))return;
		}

		for (int j = 1; j <= 3; ++j) {
			{
				uint64_t x = rnd();
				for (int k = 0; k < j; ++k)x |= rnd();
				uint64_t y = rnd() & (~x);
				assert((x & y) == 0);
				if (!onetest(x, y))return;
				if (!onetest(x, 0))return;
				if (!onetest(x, ~x))return;
			}
			{
				uint64_t x = rnd();
				for (int k = 0; k < j; ++k)x &= rnd();
				uint64_t y = rnd() & (~x);
				assert((x & y) == 0);
				if (!onetest(x, y))return;
				if (!onetest(x, 0))return;
				if (!onetest(x, ~x))return;
			}
		}

		if (!onetest(0, rnd()))return;
	}

	std::cout << "clear: test_find_functions" << std::endl;
}

void test_symmetry_functions(uint64_t seed, const int test_num = 1000000) {

	std::cout << "start: test_symmetry_functions" << std::endl;

	auto onetest = [](uint64_t x) {

		const uint64_t x060 = rotate_bitboard_060degree_clockwise_naive(x);
		const uint64_t x120 = rotate_bitboard_060degree_clockwise_naive(x060);
		const uint64_t x180 = rotate_bitboard_060degree_clockwise_naive(x120);
		const uint64_t x240 = rotate_bitboard_060degree_clockwise_naive(x180);
		const uint64_t x300 = rotate_bitboard_060degree_clockwise_naive(x240);
		const uint64_t x360 = rotate_bitboard_060degree_clockwise_naive(x300);

		if (x360 != x) {
			std::cout << "error: rotate_bitboard_060degree_clockwise_naive" << std::endl;
			std::cout << "x000 = " << std::bitset<64>(x) << std::endl;
			std::cout << "x060 = " << std::bitset<64>(x060) << std::endl;
			std::cout << "x120 = " << std::bitset<64>(x120) << std::endl;
			std::cout << "x180 = " << std::bitset<64>(x180) << std::endl;
			std::cout << "x240 = " << std::bitset<64>(x240) << std::endl;
			std::cout << "x300 = " << std::bitset<64>(x300) << std::endl;
			std::cout << "x360 = " << std::bitset<64>(x360) << std::endl;
			return false;
		}

		const uint64_t y1 = mirror_bitboard_naive(x);
		const uint64_t y2 = mirror_bitboard_naive(y1);
		if (y2 != x) {
			std::cout << "error: mirror_bitboard_naive" << std::endl;
			std::cout << "y0 = " << std::bitset<64>(x) << std::endl;
			std::cout << "y1 = " << std::bitset<64>(y1) << std::endl;
			std::cout << "y2 = " << std::bitset<64>(y2) << std::endl;
			return false;
		}

		const uint64_t x060a = rotate_bitboard_060degree_clockwise(x);
		const uint64_t x120a = rotate_bitboard_120degree_clockwise(x);
		const uint64_t x180a = rotate_bitboard_180degree_clockwise(x);
		const uint64_t x240a = rotate_bitboard_240degree_clockwise(x);
		const uint64_t x300a = rotate_bitboard_300degree_clockwise(x);

		if (x060a != x060) {
			std::cout << "error: rotate_bitboard_060degree_clockwise" << std::endl;
			std::cout << "x000  = " << std::bitset<64>(x) << std::endl;
			std::cout << "x060  = " << std::bitset<64>(x060) << std::endl;
			std::cout << "x060a = " << std::bitset<64>(x060a) << std::endl;
			return false;
		}
		if (x120a != x120) {
			std::cout << "error: rotate_bitboard_120degree_clockwise" << std::endl;
			std::cout << "x000  = " << std::bitset<64>(x) << std::endl;
			std::cout << "x120  = " << std::bitset<64>(x120) << std::endl;
			std::cout << "x120a = " << std::bitset<64>(x120a) << std::endl;
			return false;
		}
		if (x180a != x180) {
			std::cout << "error: rotate_bitboard_180degree_clockwise" << std::endl;
			std::cout << "x000  = " << std::bitset<64>(x) << std::endl;
			std::cout << "x180  = " << std::bitset<64>(x180) << std::endl;
			std::cout << "x180a = " << std::bitset<64>(x180a) << std::endl;
			return false;
		}
		if (x240a != x240) {
			std::cout << "error: rotate_bitboard_240degree_clockwise" << std::endl;
			std::cout << "x000  = " << std::bitset<64>(x) << std::endl;
			std::cout << "x240  = " << std::bitset<64>(x240) << std::endl;
			std::cout << "x240a = " << std::bitset<64>(x240a) << std::endl;
			return false;
		}
		if (x300a != x300) {
			std::cout << "error: rotate_bitboard_300degree_clockwise" << std::endl;
			std::cout << "x000  = " << std::bitset<64>(x) << std::endl;
			std::cout << "x300  = " << std::bitset<64>(x300) << std::endl;
			std::cout << "x300a = " << std::bitset<64>(x300a) << std::endl;
			return false;
		}

		const uint64_t y1a = mirror_bitboard(x);
		if (y1a != y1) {
			std::cout << "error: mirror_bitboard" << std::endl;
			std::cout << "y0  = " << std::bitset<64>(x) << std::endl;
			std::cout << "y1  = " << std::bitset<64>(y1) << std::endl;
			std::cout << "y1a = " << std::bitset<64>(y1a) << std::endl;
			return false;
		}

		return true;
	};

	for (int i = 0; i < 61; ++i) {
		const uint64_t x = 1ULL << i;
		if (!onetest(x)) {
			std::cout << "i = " << i << std::endl;
			return;
		}
	}

	std::mt19937_64 rnd(seed);
	for (int i = 0; i < test_num; ++i) {

		{
			const uint64_t x = rnd() & BB_ALL;
			if (!onetest(x)) {
				std::cout << "seed = " << seed << ", i = " << i << std::endl;
				return;
			}
		}

		for (int j = 1; j <= 3; ++j) {
			{
				uint64_t x = rnd();
				for (int k = 0; k < j; ++k)x |= rnd();
				x &= BB_ALL;
				if (!onetest(x)) {
					std::cout << "seed = " << seed << ", i = " << i << std::endl;
					return;
				}
			}
			{
				uint64_t x = rnd();
				for (int k = 0; k < j; ++k)x &= rnd();
				x &= BB_ALL;
				if (!onetest(x)) {
					std::cout << "seed = " << seed << ", i = " << i << std::endl;
					return;
				}
			}
		}
	}

	std::cout << "clear: test_symmetry_functions" << std::endl;
}

void test_pext_pdep(uint64_t seed, const int test_num = 1000000) {

	std::cout << "start: test_pext_pdep" << std::endl;

	std::mt19937_64 rnd(seed);

	const auto onetest = [](uint64_t x, uint64_t y) {
		{
			auto n1 = pext64(x, y);
			auto n2 = pext_plain(x, y);
			if (n1 != n2) {
				std::cout << "failed: x = " << x << ", y = " << y << std::endl;
				return false;
			}
		}
		{
			auto n1 = pdep64(x, y);
			auto n2 = pdep_plain(x, y);
			if (n1 != n2) {
				std::cout << "failed: x = " << x << ", y = " << y << std::endl;
				return false;
			}
		}
		return true;
	};

	for (int i = 0; i < test_num; ++i) {
		if (!onetest(rnd(), rnd()))return;
	}
	std::vector<uint64_t>s;
	for (int i = 0; i < 64; ++i) {
		uint64_t x = 1ULL << i;
		s.push_back(x);
		s.push_back(x - 1);
		s.push_back(bitwise_reverse(x));
		s.push_back(bitwise_reverse(x - 1));
	}
	for (int i = 0; i < s.size(); ++i)for (int j = 0; j < s.size(); ++j) {
		if (!onetest(s[i], s[j]))return;
		if (!onetest(s[i], rnd()))return;
		if (!onetest(rnd(), s[j]))return;
	}

	std::cout << "clear: test_pext_pdep" << std::endl;
}

void test_speed() {
	std::cout << "start: 1M random playouts from the initial position" << std::endl;
	const auto start = std::chrono::system_clock::now();
	int x[3] = {};
	for (int i = 0; i < 1000000; ++i) {
		x[std::clamp(random_play_root(i), -1, 1) + 1]++;
	}
	const auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start).count();
	std::cout << "elapsed time: " << elapsed_time << " ms" << std::endl;
	std::cout << "black win(s): " << x[2] << std::endl;
	std::cout << "white win(s): " << x[0] << std::endl;
	std::cout << "draw(s)     : " << x[1] << std::endl;

}

int main(int argc, char **argv) {

	init_rot_pos_arrays(rot060_pos, rot120_pos, rot240_pos, rot300_pos, mirror_pos);

	init_p_arrays(p_rot060, p_rot120, p_rot240, p_rot300, p_mirror, p_rot060_pop, p_rot120_pop, p_rot240_pop, p_rot300_pop, p_mirror_pop);

	if (argc >= 2) {

		test_speed();
		test_find_functions(1234);
		test_symmetry_functions(1234);
		test_search_functions(12345);
		test_pext_pdep(5555);

		test_mcts_strength(1000, 1000);

		return 0;
	}

	CLI_play_game(1000);

	return 0;
}
