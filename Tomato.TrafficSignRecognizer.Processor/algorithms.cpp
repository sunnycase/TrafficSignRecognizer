﻿//
// Traffic Sign Recognizer
// 算法
// 作者：SunnyCase
// 创建日期：2016-02-26
#include "pch.h"
#include "algorithms.h"
#include <amp_math.h>
#include <opencv2/opencv.hpp>

using namespace concurrency;
using namespace concurrency::graphics;

#define F 1000.f

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
		const auto x1 = fast_math::sqrt(x²);
		const auto x2 = -x1;

		const auto coordY = float(index[0] + y) / (float)image.extent[0];
		for (float x = x2; x <= x1; x += 0.5f)
		{
			const auto coordX = float(index[1] + x) / (float)image.extent[1];
			float theGray = Grayscale(image.sample(float_2(coordX, coordY)));
			sum += theGray;
			if (fast_math::fabs(theGray - gray) <= threshold)
				same++;
			total++;
		}
	}
	return unorm(same < (total / 2) ? 1.f : 0.f);
}

float CalculateTangent(const concurrency::graphics::texture_view<const concurrency::graphics::unorm_4, 2>& image, concurrency::index<2> index) restrict(amp)
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
	return Gy / Gx;
}

bool FitEllipse(concurrency::index<2> p1, concurrency::index<2> p2, float p1Tan, float p2Tan, const concurrency::graphics::texture_view<const concurrency::graphics::unorm, 2>& edgeView, const concurrency::graphics::texture_view<const float, 2>& tangentView, concurrency::array<uint32_t, 1>& fitsCount,
	concurrency::array<EllipseParam, 1>& ellipses, index<2>& p3) restrict(amp)
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
		// 查找次数
		const auto times = graphics::uint(MGLen / 0.5f);
		const auto step = MG / (float)times;

		const auto P1P2 = pP2 - pP1;
		auto p1p2Arctan = fast_math::atan(P1P2.y / P1P2.x);

		for (graphics::uint i = 1; i < (times - 1); i++)
		{
			const auto cntP3 = pM + (step * i);
			const auto P3Coord = coord(float_2(cntP3.x, height - cntP3.y), edgeView.extent);
			// 如果 T 点在图外则无法判断
			if (P3Coord.x < 0 || P3Coord.y < 0 || P3Coord.x > edgeView.extent[1] || P3Coord.y > edgeView.extent[0])
				continue;
			if (edgeView.sample<filter_point>(P3Coord) != 0.f)
			{
				// 判断 p3 切线是否与 P1P2平行
				auto p3Tan = tangentView.sample<filter_point>(P3Coord);
				auto m = fast_math::fabs(p1p2Arctan - fast_math::atan(p3Tan));
				if (fast_math::fabs(p1p2Arctan - fast_math::atan(p3Tan)) <= (10.f * 2 * 3.14f / 360.f))
				{
					// 平行
					pP3 = cntP3;
					p3 = index<2>(height - pP3.y, pP3.x);
					findP3 = true;
					break;
				}
			}
		}
	}
	if (findP3)
	{
		// 取 3 点为中心边长为 3 的 3 个正方形，共 27 点
		const float_2 points[] = {
			pP1 + float_2(-1, -1), pP1 + float_2(0, -1), pP1 + float_2(1, -1),
			pP1 + float_2(-1, 0),  pP1,					 pP1 + float_2(1, 0),
			pP1 + float_2(-1, 1),  pP1 + float_2(0, 1),  pP1 + float_2(1, 1),

			pP2 + float_2(-1, -1), pP2 + float_2(0, -1), pP2 + float_2(1, -1),
			pP2 + float_2(-1, 0),  pP2,					 pP2 + float_2(1, 0),
			pP2 + float_2(-1, 1),  pP2 + float_2(0, 1),  pP2 + float_2(1, 1),

			pP3 + float_2(-1, -1), pP3 + float_2(0, -1), pP3 + float_2(1, -1),
			pP3 + float_2(-1, 0),  pP3,					 pP3 + float_2(1, 0),
			pP3 + float_2(-1, 1),  pP3 + float_2(0, 1),  pP3 + float_2(1, 1)
		};
		static const graphics::uint maxPoints = sizeof(points) / sizeof(float_2);

		uint32_t rows = 0;
		float U[maxPoints][5];
		float UT[5][maxPoints];
		for (uint32_t i = 0; i < maxPoints; i++)
		{
			const auto point = points[i];
			if (edgeView.sample<filter_point>(coord(float_2(point.x, height - point.y), edgeView.extent)) != 0.f)
			{
				// x² xy y² x y
				U[rows][0] = point.x * point.x;
				U[rows][1] = point.x * point.y;
				U[rows][2] = point.y * point.y;
				U[rows][3] = point.x;
				U[rows][4] = point.y;

				UT[0][rows] = U[rows][0];
				UT[1][rows] = U[rows][1];
				UT[2][rows] = U[rows][2];
				UT[3][rows] = U[rows][3];
				UT[4][rows] = U[rows][4];
				rows++;
			}
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
				sum += UT[i][k] * F;
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
				ellipse.a = fast_math::sqrt((fast_math::pow((ellipse.E - ellipse.B * ellipse.D / (2 * ellipse.A)), 2)
					/ (4 * ellipse.C - ellipse.B * ellipse.B / ellipse.A) - F + ellipse.D *ellipse.D / (4 * ellipse.A)) / ellipse.A);
				ellipse.b = fast_math::sqrt((fast_math::pow((ellipse.E - ellipse.B * ellipse.D / (2 * ellipse.A)), 2)
					/ (4 * ellipse.C - ellipse.B * ellipse.B / ellipse.A) - F + ellipse.D *ellipse.D / (4 * ellipse.A))
					/ (ellipse.C - ellipse.B * ellipse.B / (4 * ellipse.A)));

				auto id = atomic_fetch_add(&fitsCount(0), 1);
				ellipses(id) = ellipse;
				return true;
			}
		}
	}

	return false;
}

concurrency::graphics::float_2 coord(concurrency::graphics::float_2 point, const concurrency::extent<2>& extent) restrict(cpu, amp)
{
	return float_2(point.x / extent[1], point.y / extent[0]);
}

bool Solve(float(&mat)[5][6]) restrict(cpu, amp)
{
	const auto& a = mat[0];
	for (int i = 0; i < 5; i++)
	{
		bool find = false;
		// 查找当前元的非零行
		for (int j = i; j < 5; j++)
			if (mat[j][i] != 0)
			{
				const auto factor = mat[j][i];
				for (int k = 0; k < 6; k++)
					mat[j][k] /= factor;
				swap(mat[j], mat[i]);
				find = true;
				break;
			}
		// 都是 0 则无解
		if (!find) return false;
		// 消去其他当前元不为0的行
		for (int j = 0; j < 5; j++)
		{
			if (j != i && mat[j][i] != 0)
			{
				// 要减去的factor
				const auto minusFactor = mat[j][i];
				for (int k = 0; k < 6; k++)
					mat[j][k] -= mat[i][k] * minusFactor;
			}
		}
	}
	return true;
}

bool OnEllipse(const EllipseParam & ellipse, concurrency::graphics::float_2 point, float threhold) restrict(amp)
{
	// x 方向交点 y = y0
	float x[2];
	if (Solve(ellipse.A, ellipse.B * point.y + ellipse.D, ellipse.C * point.y * point.y + ellipse.E * point.y + F, x))
	{
		if (fast_math::fabs(point.x - x[0]) <= threhold) return true;
		if (fast_math::fabs(point.x - x[1]) <= threhold) return true;
	}
	// y 方向交点 x = x0
	float y[2];
	if (Solve(ellipse.C, ellipse.B * point.x + ellipse.E, ellipse.A * point.x * point.x + ellipse.D * point.x + F, y))
	{
		if (fast_math::fabs(point.y - y[0]) <= threhold) return true;
		if (fast_math::fabs(point.y - y[1]) <= threhold) return true;
	}
	return false;
}

bool Solve(float a, float b, float c, float(&x)[2]) restrict(cpu, amp)
{
	auto tmp = b * b - 4.f * a * c;
	if (tmp < 0)
		return false;
	tmp = fast_math::sqrt(tmp);
	x[0] = (-b + tmp) / 2.f * a;
	x[1] = (-b - tmp) / 2.f * a;
	return true;
}

void HoughCircles(concurrency::array<concurrency::graphics::uint, 2>& image)
{
	std::vector<::uint> source(image.extent.size());
	copy(image, source.data());
	std::vector<byte> gray(source.size());
	std::transform(source.begin(), source.end(), gray.begin(), [](::uint color) { return (byte)(color & 0xFF); });
	cv::Mat grayMat(image.extent[0], image.extent[1], CV_8U, gray.data());

	//std::transform(grayMat.datastart, grayMat.dataend, source.begin(), [](byte color) { return ::uint(0xFF000000 | (color << 16) | (color << 8) | color);});
	//copy(source.begin(), source.end(), image);

	std::vector<std::vector<cv::Point>> contours;
	cv::findContours(grayMat, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE);
	auto p = contours.size();
}

END_NS_TSR_PRCSR