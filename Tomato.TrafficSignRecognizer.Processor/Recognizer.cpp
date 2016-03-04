//
// Traffic Sign Recognizer
// Recognizer
// 作者：SunnyCase
// 创建日期：2016-02-25
#include "pch.h"
#include "Recognizer.h"
#include <amp_math.h>
#include "algorithms.h"
#include "amprng/amp_sobol_rng.h"

using namespace NS_TSR_PRCSR;
using namespace Platform;
using namespace concurrency;
using namespace concurrency::graphics;
using namespace Windows::Foundation;
using namespace Windows::Graphics::Imaging;
using namespace WRL;

Recognizer::Recognizer(unsigned int width, unsigned int height)
	:_targetImageExtent(height, width),
	_targetImage(_targetImageExtent, 16U), _outputTex(_targetImageExtent),
	_edgeTex(_targetImageExtent, 8), _tangentTex(_targetImageExtent, 32)
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
		parallel_for_each(src.extent, [&src, targetImage](index<2> index) restrict(amp)
		{
			auto pixel = src[index];

			auto r = ((pixel & 0xFF0000) >> 16) / 255.f;
			auto g = ((pixel & 0xFF00) >> 8) / 255.f;
			auto b = (pixel & 0xFF) / 255.f;
			targetImage.set(index, unorm_4(r, g, b, 1.f));
		});
	});
}

concurrency::task<void> Recognizer::FindContours()
{
	FindEdgesAndTangent();
	texture_view<const unorm, 2> edges(_edgeTex);
	texture_view<const float, 2> tangents(_tangentTex);
	array_view<uint, 2> outputTex(_outputTex);

	return FindEllipses().then([=]
	{
		parallel_for_each(_targetImageExtent, [edges, tangents, outputTex](index<2> index) restrict(amp)
		{
			auto gray = fast_math::sin(fast_math::atan(tangents[index]));
			if (gray < 0)
			{
				auto lastGray = uint(-gray * 255);
				outputTex[index] = uint(0xFF000000 | (lastGray << 16) | (lastGray << 8));
			}
			else
			{
				auto lastGray = uint(gray * 255);
				outputTex[index] = uint(0xFF000000 | (lastGray << 8) | lastGray);
			}
		});
	});
}

void Recognizer::FindEdgesAndTangent()
{
	uint radius = 3;
	texture_view<unorm, 2> edgeView(_edgeTex);
	texture_view<float, 2> tangentView(_tangentTex);
	texture_view<const unorm_4, 2> inputTex(_targetImage);
	parallel_for_each(inputTex.extent, [inputTex, edgeView, tangentView, radius](index<2> index) restrict(amp)
	{
		auto edge = SusanTest(inputTex, index, radius, 0.1f);
		edgeView.set(index, edge);

		if (edge > 0.f)
			tangentView.set(index, CalculateTangent(inputTex, index));
	});
}

task<void> Recognizer::FindEllipses()
{
	array<uint, 2> edgePointsCount(1, 1, accelerator().default_view, access_type_read);
	array<index<2>, 2> edgePositions(_targetImageExtent);
	texture_view<const unorm, 2> edgeView(_edgeTex);

	parallel_for_each(edgeView.extent, [&edgePointsCount, &edgePositions, edgeView](index<2> index) restrict(amp)
	{
		if (edgeView[index] != 0.f)
		{
			auto idx = atomic_fetch_add(&edgePointsCount(0, 0), 1);
			auto x = idx % edgeView.extent[1];
			auto y = idx / edgeView.extent[1];
			edgePositions(y, x) = index;
		}
	});

	const auto edgeSize = extent<1>(edgePointsCount(0, 0));

	typedef sobol_rng<1>::sobol_number<float> sobol_float_number;
	sobol_rng_collection<sobol_rng<1>, 1> sc_rng(edgeSize, 1);
	concurrency::graphics::texture_view<const float, 2> tangentView(_tangentTex);
	array<float, 1> ids(edgeSize, accelerator().default_view, access_type_read);
	parallel_for_each(edgeSize / 2, [=, &edgePositions, &ids](index<1> index) restrict(amp)
	{
		const auto id1 = index[0];
		const auto x1 = id1 % edgeView.extent[1];
		const auto y1 = id1 / edgeView.extent[1];
		const auto p1Index = edgePositions(y1, x1);

		// Get the sobol RNG 
		auto rng = sc_rng[index];
		// Skip ahead to the right position
		rng.skip(sc_rng.direction_numbers(), index[0]);

		const auto id2 = int(rng.get_single(1) * edgeSize[0]);
		const auto x2 = id2 % edgeView.extent[1];
		const auto y2 = id2 / edgeView.extent[1];
		const auto p2Index = edgePositions(y2, x2);
		FitEllipse(p1Index, p2Index, tangentView[p1Index], tangentView[p2Index], edgeView, tangentView, ids[index]);
	});

	auto p = ids.data();
	auto c = std::count_if(p, p + ids.extent.size(), [](float v) { return v != 0.f;});

	float mat[5][6] = 
	{
		{ 2, 0, 0, 0, 0, 100 },
		{ 4, 1, 0, 0, 0, 100 },
		{ 0, 0, 1, 0, 0, 100 },
		{ 0, 0, 0, 1, 0, 100 },
		{ 0, 0, 0, 0, 1, 100 }
	};
	auto ret = Solve(mat);

	return task_from_result();
}
