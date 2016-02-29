//
// Traffic Sign Recognizer
// Recognizer
// 作者：SunnyCase
// 创建日期：2016-02-25
#include "pch.h"
#include "Recognizer.h"
#include <amp_math.h>
#include "algorithms.h"

using namespace NS_TSR_PRCSR;
using namespace Platform;
using namespace concurrency;
using namespace concurrency::graphics;
using namespace Windows::Foundation;
using namespace Windows::Graphics::Imaging;

Recognizer::Recognizer(unsigned int width, unsigned int height)
	:_targetImageExtent(height, width),
	_targetImage(_targetImageExtent, 16U), _outputTex(_targetImageExtent)
{

}

Windows::Foundation::IAsyncAction ^ Recognizer::SetTarget(BitmapFrame ^ frame)
{
	return create_async([=]
	{
		return PrepareTargetImage(frame);
	});
}

IAsyncOperation<bool>^ Recognizer::Recognize(Windows::Storage::Streams::IRandomAccessStream^ outputStream)
{
	return create_async([=]() -> task<bool>
	{
		return FindContours().then([=]
		{
			return create_task(BitmapEncoder::CreateAsync(BitmapEncoder::JpegEncoderId, outputStream))
				.then([=](BitmapEncoder^ encoder)
			{
				auto data = std::make_shared<std::vector<uint>>(_targetImageExtent.size());
				return copy_async(_outputTex, data->begin()).to_task().then([=]
				{
					encoder->SetPixelData(BitmapPixelFormat::Bgra8, BitmapAlphaMode::Premultiplied, _targetImageExtent[1], _targetImageExtent[0],
						96., 96., ArrayReference<byte>((byte*)data->data(), data->size() * 4));
					return create_task(encoder->FlushAsync());
				});
			});
		}).then([]
		{
			return false;
		});
	});
}

concurrency::task<void> Recognizer::PrepareTargetImage(BitmapFrame ^ frame)
{
	auto transform = ref new BitmapTransform();
	return create_task(frame->GetPixelDataAsync(BitmapPixelFormat::Bgra8, BitmapAlphaMode::Premultiplied, transform,
		ExifOrientationMode::RespectExifOrientation, ColorManagementMode::ColorManageToSRgb))
		.then([=](PixelDataProvider^ dataProvider)
	{
		auto data = dataProvider->DetachPixelData();
		array<uint, 2> src(_targetImageExtent, (uint*)data->begin(), (uint*)data->end());
		texture_view<unorm_4, 2> targetImage(_targetImage);
		parallel_for_each(src.extent, [&src, targetImage](index<2> idx) restrict(amp)
		{
			auto pixel = src[idx];

			auto r = ((pixel & 0xFF0000) >> 16) / 255.f;
			auto g = ((pixel & 0xFF00) >> 8) / 255.f;
			auto b = (pixel & 0xFF) / 255.f;
			targetImage.set(idx, unorm_4(r, g, b, 1.f));
		});
	});
}

concurrency::task<void> Recognizer::FindContours()
{
	uint radius = 3;

	int same = 0, total = 0;
	// 圆扫描线 x = n
	auto maxY = (int)radius;
	for (int y = -maxY; y <= maxY; y++)
	{
		float x² = radius * radius - y * y;
		float x1 = fast_math::sqrt(x²);
		auto x2 = -x1;

		for (float x = x2; x <= x1; x += 1.f)
			total++;
	}

	texture_view<const unorm_4, 2> inputTex(_targetImage);
	array_view<uint, 2> outputTex(_outputTex);
	parallel_for_each(inputTex.extent, [inputTex, outputTex, radius, total](index<2> idx) restrict(amp)
	{
		auto pixel = inputTex[idx];
		// 转化为灰度图
		float gray = Grayscale(pixel);

		float sum = 0.f;
		int same = 0;
		// 圆扫描线 x = n
		auto maxY = (int)radius;
		for (int y = -maxY; y <= maxY; y++)
		{
			const auto x² = float(radius * radius - y * y);
			const auto x1 = fast_math::sqrt(x²);
			const auto x2 = -x1;

			const auto coordY = float(idx[0] + y) / (float)inputTex.extent[0];
			for (float x = x2; x <= x1; x += 1.f)
			{
				const auto coordX = float(idx[1] + x) / (float)inputTex.extent[1];
				float theGray = Grayscale(inputTex.sample(float_2(coordX, coordY)));
				sum += theGray;
				auto threshold = 0.1f;
				if (fast_math::fabs(theGray - gray) <= threshold)
					same++;
			}
		}

		auto lastGray = same < (total / 2) ? 255 : 0;
		outputTex[idx] = uint(0xFF000000 | (lastGray << 16) | (lastGray << 8) | lastGray);
	});
	return task_from_result();
}
