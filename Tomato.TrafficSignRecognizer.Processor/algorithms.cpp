//
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

void HoughCircles(concurrency::array<concurrency::graphics::uint, 2>& image)
{
	std::vector<::uint> source(image.extent.size());
	copy(image, source.data());
	std::vector<byte> gray(source.size());
	std::transform(source.begin(), source.end(), gray.begin(), [](::uint color) { return (byte)(color & 0xFF);});
	cv::Mat grayMat(image.extent[0], image.extent[1], CV_8U, gray.data());

	//std::transform(grayMat.datastart, grayMat.dataend, source.begin(), [](byte color) { return ::uint(0xFF000000 | (color << 16) | (color << 8) | color);});
	//copy(source.begin(), source.end(), image);

	std::vector<std::vector<cv::Point>> contours;
	cv::findContours(grayMat, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE);
	auto p = contours.size();
}

END_NS_TSR_PRCSR