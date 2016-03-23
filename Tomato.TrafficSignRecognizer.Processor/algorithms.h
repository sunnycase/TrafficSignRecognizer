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

#define F 1000.f

concurrency::graphics::unorm Grayscale(const concurrency::graphics::unorm_4& color) restrict(cpu, amp);
concurrency::graphics::unorm Threshold(concurrency::graphics::unorm color, float threshold) restrict(cpu, amp);
concurrency::graphics::unorm SusanTest(const concurrency::graphics::texture_view<const concurrency::graphics::unorm_4, 2>& image, concurrency::index<2> index, unsigned int radius, float threshold) restrict(amp);
concurrency::graphics::unorm SusanTest(const concurrency::graphics::texture_view<const concurrency::graphics::unorm, 2>& image, concurrency::index<2> index, unsigned int radius, float threshold) restrict(amp);
concurrency::graphics::float_2 CalculateTangent(const concurrency::graphics::texture_view<const concurrency::graphics::unorm_4, 2>& image, concurrency::index<2> index) restrict(amp);
float CalculateTangent(const concurrency::graphics::texture_view<const concurrency::graphics::unorm, 2>& image, concurrency::index<2> index) restrict(amp);
bool FindEllipsePoints(concurrency::graphics::float_2(&points)[2], concurrency::graphics::float_2(&tagents)[2], const concurrency::graphics::texture_view<const concurrency::graphics::unorm, 2>& edgeView, const concurrency::graphics::texture_view<const concurrency::graphics::float_2, 2>& tangentView, EllipsePoints& ellipse) restrict(amp);
bool FitEllipse(concurrency::index<2> p1, concurrency::index<2> p2, float p1Tan, float p2Tan, const concurrency::graphics::texture_view<const concurrency::graphics::unorm, 2>& edgeView, const concurrency::graphics::texture_view<const float, 2>& tangentView, concurrency::array<uint32_t, 1>& fitsCount, concurrency::array<EllipseParam, 1>& ellipses) restrict(amp);
concurrency::graphics::float_2 coord(concurrency::graphics::float_2 point, const concurrency::extent<2>& extent) restrict(cpu, amp);
bool OnEllipse(const EllipseParam& ellipse, concurrency::graphics::float_2 point, float threhold) restrict(cpu,amp);
bool InEllipse(const EllipseParam& ellipse, concurrency::graphics::float_2 point) restrict(cpu, amp);
bool IsRed(concurrency::graphics::unorm_4 pixel) restrict(cpu, amp);
double ZernikeR(int p, int q, double rho) restrict(cpu, amp);
concurrency::graphics::double_2 ZernikeV(double r, int q, double theta) restrict(cpu, amp);
int factorial(uint32_t n) restrict(cpu, amp);
int powneg1(uint32_t n) restrict(cpu, amp);
uint32_t GetOSTUThreshold(const std::vector<uint32_t>& histGram);

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

template<uint32_t N>
bool FitEllipse(concurrency::graphics::float_2(&points)[N], double width, double height, EllipseParam& ellipse) restrict(amp)
{
	using namespace concurrency;

	double U[N][5];
	double UT[5][N];
	for (uint32_t i = 0; i < N; i++)
	{
		const auto& point = points[i];
		// x² xy y² x y
		U[i][0] = point.x * point.x;
		U[i][1] = point.x * point.y;
		U[i][2] = point.y * point.y;
		U[i][3] = point.x;
		U[i][4] = point.y;

		UT[0][i] = U[i][0];
		UT[1][i] = U[i][1];
		UT[2][i] = U[i][2];
		UT[3][i] = U[i][3];
		UT[4][i] = U[i][4];
	}

	// UT x U 形成 5 x 5 矩阵
	// UT x V 形成 5 x 1 矩阵
	double mat[5][5 + 1];
	for (uint32_t i = 0; i < 5; i++)
	{
		for (uint32_t j = 0; j < 5; j++)
		{
			double sum = 0;
			for (uint32_t k = 0; k < N; k++)
				sum += UT[i][k] * U[k][j];
			mat[i][j] = sum;
		}
		double sum = 0;
		for (uint32_t k = 0; k < N; k++)
			sum += UT[i][k] * -F;
		mat[i][5] = sum;
	}
	if (Solve(mat))
	{
		EllipseParam el{ mat[0][5], mat[1][5], mat[2][5], mat[3][5], mat[4][5] };
		// 椭圆判别式 B² - 4AC < 0
		auto value = el.B * el.B - 4.0 * el.A * el.C;
		if (value < 0.0)
		{
			el.y = (2.0 * el.A * el.E - el.B * el.D) / (el.B * el.B - 4.0 * el.A * el.C);
			el.x = (2.0 * el.C * el.D - el.B * el.E) / (el.B * el.B - 4.0 * el.A * el.C);
			if (el.x > 0.0 && el.x < width &&
				el.y > 0.0 && el.y < height)
			{
				el.a = precise_math::sqrt(2.0 * (
					(el.A * el.E * el.E - el.B * el.D * el.E + el.C * el.D * el.D) / (4.0 * el.A * el.C - el.B * el.B) - F
					) / (el.A + el.C - precise_math::sqrt(precise_math::pow(el.A - el.C, 2) + el.B * el.B)));
				el.b = precise_math::sqrt(2.0 * (
					(el.A * el.E * el.E - el.B * el.D * el.E + el.C * el.D * el.D) / (4.0 * el.A * el.C - el.B * el.B) - F
					) / (el.A + el.C + precise_math::sqrt(precise_math::pow(el.A - el.C, 2) + el.B * el.B)));

				if (el.a > 10.0 && el.b > 10.f)
				{
					el.theta = precise_math::fabs(precise_math::atan(el.B / (el.A - el.C)) / 2.0);
					el.area = uint32_t(precise_math::round(el.a * el.b * 3.1415));
					const auto major = precise_math::fmax(el.a, el.b);
					const auto minor = precise_math::fmin(el.a, el.b);
					el.length = uint32_t(precise_math::round(2.0 * 3.14 * minor + 4.0 * (major - minor)));
					el.rank = 0;

					ellipse = el;
					return true;
				}
			}
		}
	}
	return false;
}

END_NS_TSR_PRCSR