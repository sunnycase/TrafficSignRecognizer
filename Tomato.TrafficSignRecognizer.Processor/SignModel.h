//
// Traffic Sign Recognizer
// 模型
// 作者：SunnyCase
// 创建日期：2016-05-5
#include "common.h"
#include "../libsvm/svm.h"
#pragma once

DEFINE_NS_TSR_PRCSR

public ref class SignModel sealed
{
public:
	virtual ~SignModel();
internal:
	SignModel(svm_model* model);
private:
	svm_model* _model;
};

END_NS_TSR_PRCSR