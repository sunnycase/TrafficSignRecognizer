//
// Traffic Sign Recognizer
// 模型
// 作者：SunnyCase
// 创建日期：2016-05-5
#include "pch.h"
#include "SignModel.h"

using namespace NS_TSR_PRCSR;
using namespace Platform;
using namespace Platform::Collections;
using namespace concurrency;
using namespace concurrency::graphics;
using namespace Windows::Foundation;
using namespace Windows::Foundation::Collections;
using namespace Windows::Graphics::Imaging;
using namespace WRL;

SignModel::SignModel(svm_model * model)
	:_model(model)
{

}

SignModel::~SignModel()
{
	svm_free_and_destroy_model(&_model);
}
