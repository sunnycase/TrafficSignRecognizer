//
// Traffic Sign Recognizer
// 算法
// 作者：SunnyCase
// 创建日期：2016-02-26
#include "pch.h"
#include "algorithms.h"
#include "FeatureExtractor.h"
#include <amp_math.h>

using namespace concurrency;
using namespace concurrency::graphics;
using namespace Windows::Graphics::Imaging;

DEFINE_NS_TSR_PRCSR

unorm Grayscale(const unorm_4 & color) restrict(cpu, amp)
{
	return unorm(color.r * 0.299f + color.g * 0.587f + color.b * 0.114f);
}

unorm Threshold(unorm color, float threshold) restrict(cpu, amp)
{
	return unorm(color < threshold ? 0.f : 1.f);
}

concurrency::graphics::unorm SusanTest(const texture_view<const unorm_4, 2>& image, index<2> index, unsigned int radius, float threshold) restrict(amp)
{
	auto pixel = image[index];
	// 转化为灰度图
	float gray = Grayscale(pixel);

	float sum = 0.f;
	int same = 0, total = 0;
	// 圆扫描线 x = n
	auto maxY = (int)radius;
	for (int y = -maxY; y <= maxY; y++)
	{
		const auto x² = float(radius * radius - y * y);
		const auto x1 = precise_math::sqrt(x²);
		const auto x2 = -x1;

		const auto coordY = float(index[0] + y) / (float)image.extent[0];
		for (float x = x2; x <= x1; x += 0.5f)
		{
			const auto coordX = float(index[1] + x) / (float)image.extent[1];
			float theGray = Grayscale(image.sample(float_2(coordX, coordY)));
			sum += theGray;
			if (precise_math::fabs(theGray - gray) <= threshold)
				same++;
			total++;
		}
	}
	return unorm(same < (total / 2) ? 1.f : 0.f);
}

concurrency::graphics::unorm SusanTest(const texture_view<const unorm, 2>& image, index<2> index, unsigned int radius, float threshold) restrict(amp)
{
	// 获取灰度值
	float gray = image[index];

	float sum = 0.f;
	int same = 0, total = 0;
	// 圆扫描线 x = n
	auto maxY = (int)radius;
	for (int y = -maxY; y <= maxY; y++)
	{
		const auto x² = float(radius * radius - y * y);
		const auto x1 = precise_math::sqrt(x²);
		const auto x2 = -x1;

		const auto coordY = float(index[0] + y) / (float)image.extent[0];
		for (float x = x2; x <= x1; x += 0.5f)
		{
			const auto coordX = float(index[1] + x) / (float)image.extent[1];
			float theGray = image.sample(float_2(coordX, coordY));
			sum += theGray;
			if (precise_math::fabs(theGray - gray) <= threshold)
				same++;
			total++;
		}
	}
	return unorm(same < (total / 2) ? 1.f : 0.f);
}

float_2 CalculateTangent(const concurrency::graphics::texture_view<const concurrency::graphics::unorm_4, 2>& image, concurrency::index<2> index) restrict(amp)
{
	float a[9];
	a[0] = Grayscale(image.sample(float_2(float(index[1] - 1) / (float)image.extent[1], float(index[0] - 1) / (float)image.extent[0])));
	a[1] = Grayscale(image.sample(float_2(float(index[1]) / (float)image.extent[1], float(index[0] - 1) / (float)image.extent[0])));
	a[2] = Grayscale(image.sample(float_2(float(index[1] + 1) / (float)image.extent[1], float(index[0] - 1) / (float)image.extent[0])));

	a[3] = Grayscale(image.sample(float_2(float(index[1] - 1) / (float)image.extent[1], float(index[0]) / (float)image.extent[0])));
	a[5] = Grayscale(image.sample(float_2(float(index[1] + 1) / (float)image.extent[1], float(index[0]) / (float)image.extent[0])));

	a[6] = Grayscale(image.sample(float_2(float(index[1] - 1) / (float)image.extent[1], float(index[0] + 1) / (float)image.extent[0])));
	a[7] = Grayscale(image.sample(float_2(float(index[1]) / (float)image.extent[1], float(index[0] + 1) / (float)image.extent[0])));
	a[8] = Grayscale(image.sample(float_2(float(index[1] + 1) / (float)image.extent[1], float(index[0] + 1) / (float)image.extent[0])));

	auto Gx = (a[6] + 2 * a[7] + a[8]) - (a[0] + 2 * a[1] + a[2]);
	auto Gy = (a[2] + 2 * a[5] + a[8]) - (a[0] + 2 * a[3] + a[6]);
	return float_2(Gx, Gy);
}

float CalculateTangent(const concurrency::graphics::texture_view<const concurrency::graphics::unorm, 2>& image, concurrency::index<2> index) restrict(amp)
{
	float a[9];
	a[0] = image.sample(float_2(float(index[1] - 1) / (float)image.extent[1], float(index[0] - 1) / (float)image.extent[0]));
	a[1] = image.sample(float_2(float(index[1]) / (float)image.extent[1], float(index[0] - 1) / (float)image.extent[0]));
	a[2] = image.sample(float_2(float(index[1] + 1) / (float)image.extent[1], float(index[0] - 1) / (float)image.extent[0]));

	a[3] = image.sample(float_2(float(index[1] - 1) / (float)image.extent[1], float(index[0]) / (float)image.extent[0]));
	a[5] = image.sample(float_2(float(index[1] + 1) / (float)image.extent[1], float(index[0]) / (float)image.extent[0]));

	a[6] = image.sample(float_2(float(index[1] - 1) / (float)image.extent[1], float(index[0] + 1) / (float)image.extent[0]));
	a[7] = image.sample(float_2(float(index[1]) / (float)image.extent[1], float(index[0] + 1) / (float)image.extent[0]));
	a[8] = image.sample(float_2(float(index[1] + 1) / (float)image.extent[1], float(index[0] + 1) / (float)image.extent[0]));

	auto Gx = (a[6] + 2 * a[7] + a[8]) - (a[0] + 2 * a[1] + a[2]);
	auto Gy = (a[2] + 2 * a[5] + a[8]) - (a[0] + 2 * a[3] + a[6]);
	return Gy / Gx;
}

float_2 FixPoint(float_2 point) restrict(amp)
{
	return float_2(precise_math::round(point.x), precise_math::round(point.y));
}

bool FindEllipsePoints(float_2(&points)[2], float_2(&tagents)[2], const concurrency::graphics::texture_view<const concurrency::graphics::unorm, 2>& edgeView, const concurrency::graphics::texture_view<const float_2, 2>& tangentView, EllipsePoints& ellipse) restrict(amp)
{
	// 如果两切线平行则无法判断
	if (tagents[0] == tagents[1]) return false;

	const float height = edgeView.extent[0];
	const auto p1Tan = tagents[0].y / tagents[0].x;
	const auto p2Tan = tagents[1].y / tagents[1].x;

	const auto pP1 = points[0];
	const auto pP2 = points[1];
	// 求交点(椭圆的极) T
	const auto b1 = pP1.y - pP1.x * p1Tan;
	const auto b2 = pP2.y - pP1.x * p2Tan;
	const auto pTx = (b1 - b2) / (p2Tan - p1Tan);
	const auto pT = float_2(pTx, pTx * p1Tan + b1);
	if (pT.x < 0.f || pT.x > edgeView.extent[1] || pT.y < 0.f || pT.y > edgeView.extent[0])
		return false;
	// P1P2 的中点 M
	const auto pM = (pP1 + pP2) / 2.f;
	// MT 的中点 G
	const auto pG = (pM + pT) / 2.f;

	bool findP3 = false;
	float_2 pP3;
	// 线段 MG 上搜索 P3
	{
		const auto MG = pG - pM;
		const auto MGLen = fast_math::sqrt(MG.x * MG.x + MG.y * MG.y);
		if (MGLen < 5.f || MGLen > 50.f) return false;
		// 查找次数
		const auto times = int(MGLen / 1.f);
		const auto step = MG / (float)times;

		const auto P1P2 = pP2 - pP1;
		auto p1p2Arctan = fast_math::atan2(P1P2.y, P1P2.x);

		for (int i = 0; i < times; i++)
		{
			const auto cntP3 = FixPoint(pM + (step * i));
			const auto P3Coord = coord(float_2(cntP3.x + 0.5f, height - cntP3.y + 0.5f), edgeView.extent);
			// 如果 T 点在图外则无法判断
			if (P3Coord.x < 0 || P3Coord.y < 0 || P3Coord.x > 1.f || P3Coord.y > 1.f)
				continue;
			if (edgeView.sample<filter_point>(P3Coord) > 0.f)
			{
				// 判断 p3 切线是否与 P1P2平行
				auto p3Tan = tangentView.sample<filter_point>(P3Coord);
				const auto threhold = (5 * 2 * 3.14f / 360.f);
				if (fast_math::fabs(p1p2Arctan - precise_math::atan2(p3Tan.y, p3Tan.x)) <= threhold
					|| fast_math::fabs(3.14f - p1p2Arctan + precise_math::atan2(p3Tan.y, p3Tan.x)) <= threhold)
				{
					// 平行
					pP3 = cntP3;
					ellipse = { float_2(pP1.x, pP1.y), float_2(pP2.x, pP2.y), float_2(pP3.x, pP3.y)};
					return true;
				}
			}
		}
	}
	return false;
}

bool FindPoint(float_2 p1, float_2 p2, const texture_view<const unorm, 2>& edgeView, float_2& point) restrict(amp)
{
	const float height = edgeView.extent[0];
	const auto P1P2 = p2 - p1;
	const auto P1P2Len = precise_math::sqrt(P1P2.x * P1P2.x + P1P2.y * P1P2.y);
	if (P1P2Len < 1.f) return false;
	// 查找次数
	const auto times = int(P1P2Len / 1.f);
	const auto step = P1P2 / (float)times;

	for (int i = 0; i < times; i++)
	{
		const auto cntPoint = FixPoint(p1 + (step * i));
		const auto P3Coord = coord(float_2(cntPoint.x, height - cntPoint.y), edgeView.extent);
		if (P3Coord.x < 0 || P3Coord.y < 0 || P3Coord.x > 1.f || P3Coord.y > 1.f)
			continue;
		if (edgeView.sample<filter_point>(P3Coord) == 1.f)
		{
			point = cntPoint;
			return true;
		}
	}
	return false;
}

template<uint32_t N>
void FindPoint(float_2(&points)[N], uint32_t& offset, float_2 p1, float_2 pM, float_2 p3, const texture_view<const unorm, 2>& edgeView) restrict(amp)
{
	const auto P1M = pM - p1;
	const auto P1P3 = p3 - p1;

	const auto P1MLen = precise_math::sqrt(P1M.x * P1M.x + P1M.y * P1M.y);
	if (P1MLen < 1.f) return;
	const auto P1P3Len = precise_math::sqrt(P1P3.x * P1P3.x + P1P3.y * P1P3.y);
	if (P1P3Len < 1.f) return;
	const auto times = int(precise_math::fmin(P1MLen, P1P3Len) / 1.5f);
	const auto step1 = P1M / (float)times;
	const auto step2 = P1P3 / (float)times;

	for (int i = 0; i < times && offset < N; i++)
	{
		const auto cntP1M = FixPoint(p1 + (step1 * i));
		const auto cntP1P3 = FixPoint(p1 + (step2 * i));
		const auto cntQ = 2.f * cntP1P3 - cntP1M;
		float_2 point;
		if (FindPoint(cntP1P3, cntQ, edgeView, point))
		{
			const float_2 nearPoints[] =
			{
				point + float_2(-1, -1), point + float_2(+0, -1), point + float_2(+1, -1),
				point + float_2(-1, +0), point + float_2(+0, +0), point + float_2(+1, +0),
				point + float_2(-1, +1), point + float_2(+0, +1), point + float_2(+1, +1),
			};
			for (uint32_t k = 0; k < 9; k++)
			{
				point = nearPoints[k];
				const auto pCoord = coord(float_2(point.x, edgeView.extent[0] - point.y), edgeView.extent);
				if (edgeView.sample<filter_point>(pCoord) == 1.f)
				{
					bool alreadyHas = false;
					for (uint32_t j = 0; j < offset; j++)
					{
						if (points[j] == point)
						{
							alreadyHas = true;
							break;
						}
					}
					if (!alreadyHas)
						points[offset++] = point;
				}
			}
		}
	}
}

bool FitEllipse(concurrency::index<2> p1, concurrency::index<2> p2, float p1Tan, float p2Tan, const concurrency::graphics::texture_view<const concurrency::graphics::unorm, 2>& edgeView, const concurrency::graphics::texture_view<const float, 2>& tangentView, concurrency::array<uint32_t, 1>& fitsCount,
	concurrency::array<EllipseParam, 1>& ellipses) restrict(amp)
{
	// 如果两切线平行则无法判断
	if (p1Tan == p2Tan) return false;

	const float height = edgeView.extent[0];

	const auto pP1 = float_2(p1[1], height - p1[0]);
	const auto pP2 = float_2(p2[1], height - p2[0]);
	// 求交点(椭圆的极) T
	const auto b1 = pP1.y - pP1.x * p1Tan;
	const auto b2 = pP2.y - pP1.x * p2Tan;
	const auto pTx = (b1 - b2) / (p2Tan - p1Tan);
	const auto pT = float_2(pTx, pTx * p1Tan + b1);
	if (pT.x < 0.f || pT.x > edgeView.extent[1] || pT.y < 0.f || pT.y > edgeView.extent[0])
		return false;
	// P1P2 的中点 M
	const auto pM = (pP1 + pP2) / 2.f;
	// MT 的中点 G
	const auto pG = (pM + pT) / 2.f;

	bool findP3 = false;
	float_2 pP3;
	// 线段 MG 上搜索 P3
	{
		const auto MG = pG - pM;
		const auto MGLen = precise_math::sqrt(MG.x * MG.x + MG.y * MG.y);
		if (MGLen < 1.f) return false;
		// 查找次数
		const auto times = int(MGLen / 1.f);
		const auto step = MG / (float)times;

		const auto P1P2 = pP2 - pP1;
		auto p1p2Arctan = precise_math::atan(P1P2.y / P1P2.x);

		for (int i = 0; i < times; i++)
		{
			const auto cntP3 = FixPoint(pM + (step * i));
			const auto P3Coord = coord(float_2(cntP3.x, height - cntP3.y), edgeView.extent);
			// 如果 T 点在图外则无法判断
			if (P3Coord.x < 0 || P3Coord.y < 0 || P3Coord.x > 1.f || P3Coord.y > 1.f)
				continue;
			if (edgeView.sample<filter_point>(P3Coord) == 1.f)
			{
				// 判断 p3 切线是否与 P1P2平行
				auto p3Tan = tangentView.sample<filter_point>(P3Coord);
				auto m = precise_math::fabs(p1p2Arctan - precise_math::atan(p3Tan));
				if (precise_math::fabs(p1p2Arctan - precise_math::atan(p3Tan)) <= (10.f * 2 * 3.14f / 360.f))
				{
					// 平行
					pP3 = cntP3;
					{
						// 取 3 点为中心边长为 3 的 3 个正方形，共 27 点
						static const graphics::uint maxPoints = 120;
						float_2 points[maxPoints];
						points[0] = pP1;
						points[1] = pP2;
						points[2] = pP3;
						uint32_t rows = 3;
						FindPoint(points, rows, pP1, pM, pP3, edgeView);
						FindPoint(points, rows, pP2, pM, pP3, edgeView);

						float U[maxPoints][5];
						float UT[5][maxPoints];
						for (uint32_t i = 0; i < rows; i++)
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
						float mat[5][5 + 1];
						for (uint32_t i = 0; i < 5; i++)
						{
							for (uint32_t j = 0; j < 5; j++)
							{
								float sum = 0;
								for (uint32_t k = 0; k < rows; k++)
									sum += UT[i][k] * U[k][j];
								mat[i][j] = sum;
							}
							float sum = 0;
							for (uint32_t k = 0; k < rows; k++)
								sum += UT[i][k] * -F;
							mat[i][5] = sum;
						}
						if (auto solved = Solve(mat))
						{
							EllipseParam ellipse{ mat[0][5], mat[1][5], mat[2][5], mat[3][5], mat[4][5] };
							// 椭圆判别式 B² - 4AC < 0
							auto value = ellipse.B * ellipse.B - 4.f * ellipse.A * ellipse.C;
							if (value < 0.f)
							{
								ellipse.y = height - -(ellipse.E - ellipse.B * ellipse.D / (2 * ellipse.A))
									/ (2 * ellipse.C - ellipse.B * ellipse.B / (2 * ellipse.A));
								ellipse.x = -(ellipse.B * ellipse.y + ellipse.D) / (2 * ellipse.A);
								if (ellipse.x > 0.f && ellipse.x < edgeView.extent[1] &&
									ellipse.y > 0.f && ellipse.y < edgeView.extent[0])
								{
									ellipse.a = precise_math::sqrt((precise_math::pow((ellipse.E - ellipse.B * ellipse.D / (2 * ellipse.A)), 2)
										/ (4 * ellipse.C - ellipse.B * ellipse.B / ellipse.A) - F + ellipse.D *ellipse.D / (4 * ellipse.A)) / ellipse.A);
									ellipse.b = precise_math::sqrt((precise_math::pow((ellipse.E - ellipse.B * ellipse.D / (2 * ellipse.A)), 2)
										/ (4 * ellipse.C - ellipse.B * ellipse.B / ellipse.A) - F + ellipse.D *ellipse.D / (4 * ellipse.A))
										/ (ellipse.C - ellipse.B * ellipse.B / (4 * ellipse.A)));
									//for (uint32_t i = 0; i < rows; i++)
									//	ellipse.points[i] = float_2(points[i].x, height - points[i].y);
									//ellipse.p1 = float_2(pP1.x, height - pP1.y);
									//ellipse.p2 = float_2(pP2.x, height - pP2.y);
									//ellipse.p3 = float_2(pP3.x, height - pP3.y);
									auto id = atomic_fetch_add(&fitsCount(0), 1);
									ellipses(id) = ellipse;
									break;
								}
							}
						}
					}
				}
			}
		}
	}
	return false;
}

concurrency::graphics::float_2 coord(concurrency::graphics::float_2 point, const concurrency::extent<2>& extent) restrict(cpu, amp)
{
	return float_2((point.x + 0.5f) / (float)extent[1], (point.y + 0.5f) / (float)extent[0]);
}

bool OnEllipse(const EllipseParam & ellipse, concurrency::graphics::float_2 point, float threhold) restrict(cpu, amp)
{
	// x 方向交点 y = y0
	float x[2];
	if (Solve(ellipse.A, ellipse.B * point.y + ellipse.D, ellipse.C * point.y * point.y + ellipse.E * point.y + F, x))
	{
		if (precise_math::fabs(point.x - x[0]) <= threhold) return true;
		if (precise_math::fabs(point.x - x[1]) <= threhold) return true;
	}
	// y 方向交点 x = x0
	float y[2];
	if (Solve(ellipse.C, ellipse.B * point.x + ellipse.E, ellipse.A * point.x * point.x + ellipse.D * point.x + F, y))
	{
		if (precise_math::fabs(point.y - y[0]) <= threhold) return true;
		if (precise_math::fabs(point.y - y[1]) <= threhold) return true;
	}
	return false;
}

bool InEllipse(const EllipseParam & ellipse, concurrency::graphics::float_2 point) restrict(cpu, amp)
{
	auto value = ellipse.A * point.x * point.x + ellipse.B * point.x * point.y + ellipse.C * point.y * point.y + ellipse.D * point.x +
		ellipse.E * point.y + F;
	return value <= 0.f;
}

bool IsRed(concurrency::graphics::unorm_4 pixel) restrict(cpu, amp)
{
	const auto max = precise_math::fmax(pixel.r, precise_math::fmax(pixel.g, pixel.b));
	const auto min = precise_math::fmin(pixel.r, precise_math::fmin(pixel.g, pixel.b));
	const auto V = max;
	const auto S = (max - min) / max;
	float H;
	if (pixel.r == max)
		H = (pixel.g - pixel.b) / (max - min) * 60;
	if (pixel.g == max) H = 120 + (pixel.b - pixel.r) / (max - min) * 60;
	if (pixel.b == max) H = 240 + (pixel.r - pixel.g) / (max - min) * 60;
	if (H < 0) H = H + 360;
	const auto y = Grayscale(pixel);
	return ((H >= 0 && H <= 16) || (H >= 315 && H <= 360)) && (S > 0.4) && (y > 0.1 && y < 0.8);
}

double ZernikeR(int p, int q, double rho) restrict(cpu, amp)
{
	const auto sEnd = (p - abs(q)) / 2;

	double r = 0.0;
	for (int s = 0; s <= sEnd; s++)
		r += powneg1(s) * (double)factorial(p - s) /
		float(factorial(s) * factorial((p + abs(q)) / 2 - s) * factorial((p - abs(q)) / 2 - s)) *
		precise_math::pow(rho, p - 2 * s) * 255.0;
	return r;
}

concurrency::graphics::double_2 ZernikeV(double r, int q, double theta) restrict(cpu, amp)
{
	return double_2(r * precise_math::cos(q * theta),
		r * precise_math::sin(q * theta));
}

int factorial(uint32_t n) restrict(cpu, amp)
{
	if (n <= 10)
	{
		const int factorials[11] = { 1 , 1 , 2 , 6 , 24 , 120 , 720 , 5040 , 40320 , 362880 , 39916800 };
		return factorials[n];
	}

	int m = 1;
	for (uint32_t i = 2; i <= n; i++)
	{
		m *= i;
	}
	return m;
}

int powneg1(uint32_t n) restrict(cpu, amp)
{
	return (n % 2) ? -1 : 1;
}

uint32_t GetOSTUThreshold(const std::vector<uint32_t>& HistGram)
{
	int X, Y, Amount = 0;
	int PixelBack = 0, PixelFore = 0, PixelIntegralBack = 0, PixelIntegralFore = 0, PixelIntegral = 0;
	double OmegaBack, OmegaFore, MicroBack, MicroFore, SigmaB, Sigma;              // 类间方差;
	int MinValue, MaxValue;
	int Threshold = 0;

	for (MinValue = 0; MinValue < 256 && HistGram[MinValue] == 0; MinValue++);
	for (MaxValue = 255; MaxValue > MinValue && HistGram[MinValue] == 0; MaxValue--);
	if (MaxValue == MinValue) return MaxValue;          // 图像中只有一个颜色             
	if (MinValue + 1 == MaxValue) return MinValue;      // 图像中只有二个颜色

	for (Y = MinValue; Y <= MaxValue; Y++) Amount += HistGram[Y];        //  像素总数

	PixelIntegral = 0;
	for (Y = MinValue; Y <= MaxValue; Y++) PixelIntegral += HistGram[Y] * Y;
	SigmaB = -1;
	for (Y = MinValue; Y < MaxValue; Y++)
	{
		PixelBack = PixelBack + HistGram[Y];
		PixelFore = Amount - PixelBack;
		OmegaBack = (double)PixelBack / Amount;
		OmegaFore = (double)PixelFore / Amount;
		PixelIntegralBack += HistGram[Y] * Y;
		PixelIntegralFore = PixelIntegral - PixelIntegralBack;
		MicroBack = (double)PixelIntegralBack / PixelBack;
		MicroFore = (double)PixelIntegralFore / PixelFore;
		Sigma = OmegaBack * OmegaFore * (MicroBack - MicroFore) * (MicroBack - MicroFore);
		if (Sigma > SigmaB)
		{
			SigmaB = Sigma;
			Threshold = Y;
		}
	}
	return Threshold;
}

concurrency::task<std::array<double, 11>> GetFeature(Windows::Storage::StorageFile ^ file)
{
	auto decoder = co_await create_task(BitmapDecoder::CreateAsync(co_await create_task(file->OpenReadAsync())));
	auto frame = co_await create_task(decoder->GetFrameAsync(0));
	auto ext = ref new FeatureExtractor(frame->OrientedPixelWidth, frame->OrientedPixelHeight);
	co_await create_task(ext->SetTarget(frame));
	co_await create_task(ext->Recognize());
	auto rawFeatures = co_await create_task(ext->CaculateZernikes());
	for (auto&& el : rawFeatures)
	{
		std::array<double, 11> arr{ el->GetAt(0).z, el->GetAt(1).z, el->GetAt(2).z, el->GetAt(3).z, el->GetAt(4).z, el->GetAt(5).z,
			el->GetAt(6).z, el->GetAt(7).z, el->GetAt(8).z, el->GetAt(9).z, el->GetAt(10).z };
		return arr;
	}
	return std::array<double, 11>{};
}

bool Solve(float a, float b, float c, float(&x)[2]) restrict(cpu, amp)
{
	auto tmp = b * b - 4.f * a * c;
	if (tmp < 0)
		return false;
	tmp = precise_math::sqrt(tmp);
	x[0] = (-b + tmp) / (2.f * a);
	x[1] = (-b - tmp) / (2.f * a);
	return true;
}

END_NS_TSR_PRCSR