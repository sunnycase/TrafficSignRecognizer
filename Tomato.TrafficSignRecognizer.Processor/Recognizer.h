//
// Traffic Sign Recognizer
// Recognizer
// 作者：SunnyCase
// 创建日期：2016-02-25
#include "common.h"
#include <amp.h>
#include <amp_graphics.h>
#include <ppltasks.h>
#pragma once

DEFINE_NS_TSR_PRCSR

public ref class Recognizer sealed
{
public:
	Recognizer(unsigned int width, unsigned int height);

	Windows::Foundation::IAsyncAction^ SetTarget(Windows::Graphics::Imaging::BitmapFrame ^ frame);
	Windows::Foundation::IAsyncOperation<bool>^ Recognize(Windows::Storage::Streams::IRandomAccessStream^ outputStream);
private:
	concurrency::task<void> PrepareTargetImage(Windows::Graphics::Imaging::BitmapFrame ^ frame);
	concurrency::task<void> FindContours();
	concurrency::graphics::texture<concurrency::graphics::unorm, 2> FindEdges();
private:
	concurrency::extent<2> _targetImageExtent;
	concurrency::graphics::texture<concurrency::graphics::unorm_4, 2> _targetImage;
	concurrency::array<concurrency::graphics::uint, 2> _outputTex;
};

END_NS_TSR_PRCSR