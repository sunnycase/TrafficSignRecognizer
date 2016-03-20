//
// Traffic Sign Recognizer
// 算法
// 作者：SunnyCase
// 创建日期：2016-02-26
#include "common.h"
#include <amp.h>
#include <amp_graphics.h>
#include "models.h"

DEFINE_NS_TSR_PRCSR

concurrency::graphics::unorm Grayscale(const concurrency::graphics::unorm_4& color) restrict(cpu, amp);
concurrency::graphics::unorm Threshold(concurrency::graphics::unorm color, float threshold) restrict(cpu, amp);
concurrency::graphics::unorm SusanTest(const concurrency::graphics::texture_view<const concurrency::graphics::unorm_4, 2>& image, concurrency::index<2> index, unsigned int radius, float threshold) restrict(amp);
concurrency::graphics::unorm SusanTest(const concurrency::graphics::texture_view<const concurrency::graphics::unorm, 2>& image, concurrency::index<2> index, unsigned int radius, float threshold) restrict(amp);
float CalculateTangent(const concurrency::graphics::texture_view<const concurrency::graphics::unorm_4, 2>& image, concurrency::index<2> index) restrict(amp);
float CalculateTangent(const concurrency::graphics::texture_view<const concurrency::graphics::unorm, 2>& image, concurrency::index<2> index) restrict(amp);
void FindEllipsePoints(concurrency::index<2> p1, concurrency::index<2> p2, float p1Tan, float p2Tan, const concurrency::graphics::texture_view<const concurrency::graphics::unorm, 2>& edgeView, const concurrency::graphics::texture_view<const float, 2>& tangentView, concurrency::array<uint32_t, 1>& fitsCount, concurrency::array<EllipsePoints, 1>& ellipses) restrict(amp);
bool FitEllipse(concurrency::index<2> p1, concurrency::index<2> p2, float p1Tan, float p2Tan, const concurrency::graphics::texture_view<const concurrency::graphics::unorm, 2>& edgeView, const concurrency::graphics::texture_view<const float, 2>& tangentView, concurrency::array<uint32_t, 1>& fitsCount, concurrency::array<EllipseParam, 1>& ellipses) restrict(amp);
bool FitEllipse(concurrency::graphics::float_2 (&points)[5], float width, float height, concurrency::array<uint32_t, 1>& fitsCount, concurrency::array<EllipseParam, 1>& ellipses) restrict(amp);
concurrency::graphics::float_2 coord(concurrency::graphics::float_2 point, const concurrency::extent<2>& extent) restrict(cpu, amp);
bool OnEllipse(const EllipseParam& ellipse, concurrency::graphics::float_2 point, float threhold) restrict(cpu,amp);
bool InEllipse(const EllipseParam& ellipse, concurrency::graphics::float_2 point) restrict(cpu, amp);
bool IsRed(concurrency::graphics::unorm_4 pixel) restrict(cpu, amp);
double ZernikeR(int p, int q, double rho) restrict(cpu, amp);
concurrency::graphics::double_2 ZernikeV(double r, int q, double theta) restrict(cpu, amp);
int factorial(uint32_t n) restrict(cpu, amp);
int powneg1(uint32_t n) restrict(cpu, amp);

// 解一元二次方程
bool Solve(float a, float b, float c, float(&x)[2]) restrict(cpu, amp);

template<typename T, uint32_t N>
void swap(T(&v1)[N], T(&v2)[N]) restrict(cpu, amp)
{
	T tmp[N];
	for (uint32_t i = 0; i < N; i++)
		tmp[i] = v1[i];
	for (uint32_t i = 0; i < N; i++)
		v1[i] = v2[i];
	for (uint32_t i = 0; i < N; i++)
		v2[i] = tmp[i];
}

template<typename T>
void swap(T& left, T& right) restrict(amp)
{
	T tmp = left;
	left = right;
	right = tmp;
}

template<typename T, uint32_t N>
void sort(T(&arr)[N]) restrict(cpu, amp)
{
	for (int j = 0; j < N; j++)
	{
		for (int k = 0; k < N - j; k++)
		{
			if (arr[k] > arr[k + 1])
				swap(arr[k], arr[k + 1]);
		}
	}
}

template<typename T>
T abs(T value) restrict(cpu, amp)
{
	return value < 0 ? -value : value;
}

template<typename T, uint32_t M, uint32_t N = M + 1>
bool Solve(T(&mat)[M][N]) restrict(cpu, amp)
{
	for (int i = 0; i < M; i++)
	{
		bool find = false;
		// 查找当前元的非零行
		for (int j = i; j < M; j++)
			if (mat[j][i] != 0)
			{
				const auto factor = mat[j][i];
				for (int k = 0; k < N; k++)
					mat[j][k] /= factor;
				swap(mat[j], mat[i]);
				find = true;
				break;
			}
		// 都是 0 则无解
		if (!find) return false;
		// 消去其他当前元不为0的行
		for (int j = 0; j < M; j++)
		{
			if (j != i && mat[j][i] != 0)
			{
				// 要减去的factor
				const auto minusFactor = mat[j][i];
				for (int k = 0; k < N; k++)
					mat[j][k] -= mat[i][k] * minusFactor;
			}
		}
	}
	return true;
}

END_NS_TSR_PRCSR