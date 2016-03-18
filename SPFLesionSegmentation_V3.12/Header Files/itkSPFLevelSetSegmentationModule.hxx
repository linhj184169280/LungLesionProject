/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkSPFLevelSetSegmentationModule.hxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkSPFLevelSetSegmentationModule_hxx
#define __itkSPFLevelSetSegmentationModule_hxx

#include "itkSPFLevelSetSegmentationModule.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkSigmoidImageFilter.h"

namespace itk
{

/**
 * Constructor
 */
template <unsigned int NDimension>
SPFLevelSetSegmentationModule<NDimension>
::SPFLevelSetSegmentationModule()
{
  this->SetNumberOfRequiredInputs( 2 );
  this->SetNumberOfRequiredOutputs( 1 );

  typename OutputSpatialObjectType::Pointer outputObject = OutputSpatialObjectType::New();

  this->ProcessObject::SetNthOutput( 0, outputObject.GetPointer() );

  this->m_MaximumNumberOfIterations = 180;
  this->m_MinChangeNumber = 3;

  this->m_Sigma = 0.3;
  this->m_NeiborSize = 5;
  this->m_delta = 1.0;
  this->m_Propagation= 10.0;
//  this->m_ZeroSetInputImage = NULL;
  this->m_InvertOutputIntensities = true;

}


/**
 * Destructor
 */
template <unsigned int NDimension>
SPFLevelSetSegmentationModule<NDimension>
::~SPFLevelSetSegmentationModule()
{
}


/**
 * PrintSelf
 */
template <unsigned int NDimension>
void
SPFLevelSetSegmentationModule<NDimension>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
//  os << indent << "PropagationScaling = " << this->m_PropagationScaling << std::endl;
//  os << indent << "CurvatureScaling = " << this->m_CurvatureScaling << std::endl;
//  os << indent << "AdvectionScaling = " << this->m_AdvectionScaling << std::endl;
//  os << indent << "MaximumRMSError = " << this->m_MaximumRMSError << std::endl;
//  os << indent << "MaximumNumberOfIterations = " << this->m_MaximumNumberOfIterations << std::endl;
}


/**
 * This method is intended to be used only by the subclasses to extract the
 * input image from the input SpatialObject.
 */
template <unsigned int NDimension>
const typename SPFLevelSetSegmentationModule<NDimension>::InputSpatialObjectType *
SPFLevelSetSegmentationModule<NDimension>
::GetInternalInputLandmarks() const
{
	const InputSpatialObjectType * inputObject =
		dynamic_cast< const InputSpatialObjectType * >( this->GetInput() );

	return inputObject;
}


/**
 * This method is intended to be used only by the subclasses to extract the
 * input feature image from the input feature SpatialObject.
 */
template <unsigned int NDimension>
const typename SPFLevelSetSegmentationModule<NDimension>::FeatureImageType *
SPFLevelSetSegmentationModule<NDimension>
::GetInternalFeatureImage() const
{
  const FeatureSpatialObjectType * featureObject =
    dynamic_cast< const FeatureSpatialObjectType * >( this->GetFeature() );

  const FeatureImageType * featureImage = featureObject->GetImage();

  return featureImage;
}


/**
 * This method is intended to be used only by the subclasses to insert the
 * output image as cargo of the output spatial object.
 */
template <unsigned int NDimension>
void
SPFLevelSetSegmentationModule<NDimension>
::PackOutputImageInOutputSpatialObject( OutputImageType * image )
{
  typename OutputImageType::Pointer outputImage = image;

  if( this->m_InvertOutputIntensities )
    {
    typedef MinimumMaximumImageCalculator< OutputImageType > CalculatorType;
    typename CalculatorType::Pointer calculator = CalculatorType::New();
    calculator->SetImage( outputImage );
    calculator->Compute();
    typedef IntensityWindowingImageFilter< OutputImageType, OutputImageType > RescaleFilterType;
    typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
    rescaler->SetInput( outputImage );
    rescaler->SetWindowMinimum( calculator->GetMinimum() );
    rescaler->SetWindowMaximum( calculator->GetMaximum() ); 
    rescaler->SetOutputMinimum( 4.0 ); // Note that the values must be [4:-4] here to 
    rescaler->SetOutputMaximum( -4.0 ); // make sure that we invert and not just rescale.
    rescaler->InPlaceOn();
    rescaler->Update();
    outputImage = rescaler->GetOutput();
    }

  outputImage->DisconnectPipeline();

  OutputSpatialObjectType * outputObject =
    dynamic_cast< OutputSpatialObjectType * >(this->ProcessObject::GetOutput(0));

  outputObject->SetImage( outputImage );
}


/**
 * Generate Data
 */
template <unsigned int NDimension>
void
SPFLevelSetSegmentationModule<NDimension>
::GenerateData()
{
	FeatureImageType *featureImage = const_cast<FeatureImageType *>(this->GetInternalFeatureImage());
	
	////////////////////////////////////BEGIN//////////////////////
	ImageFileWriter<OutputImageType>::Pointer writer = ImageFileWriter<OutputImageType>::New();
	writer->SetInput( featureImage );
	writer->SetFileName("E:\\result\\SPF_feature_Befor_PreProcess.mha");
	try
	{
	//	writer->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;

	}
	
	///////////////////////////////////////END////////////////////////////
	typedef DiscreteGaussianImageFilter< 
		InternalImageType, InternalImageType >		  GaussianFilterType;

	typedef SigmoidImageFilter<
		InternalImageType, OutputImageType >                 SigmoidFilterType;

	GaussianFilterType::Pointer   m_PreGaussianFilter = GaussianFilterType::New();	//预处理的高斯
	SigmoidFilterType::Pointer    m_PreSigmoidFilter = SigmoidFilterType::New();

	m_PreGaussianFilter->SetInput( featureImage );
	m_PreGaussianFilter->SetUseImageSpacingOn();
	m_PreGaussianFilter->SetMaximumKernelWidth( 5 );
	m_PreGaussianFilter->SetVariance( 0.35 );

	typedef MinimumMaximumImageCalculator< OutputImageType > CalculatorType;
	typename CalculatorType::Pointer calculator = CalculatorType::New();
	calculator->SetImage( featureImage );
	calculator->Compute();


//	m_PreSigmoidFilter->SetInput( m_PreGaussianFilter->GetOutput() );
	m_PreSigmoidFilter->SetInput( featureImage );
	m_PreSigmoidFilter->SetAlpha( 100.0 );
	m_PreSigmoidFilter->SetBeta( -250.0 );
	m_PreSigmoidFilter->SetOutputMinimum( 0.0 );
	m_PreSigmoidFilter->SetOutputMaximum( 1.0 );
//	m_PreSigmoidFilter->SetOutputMinimum( (float)calculator->GetMinimum() );
//	m_PreSigmoidFilter->SetOutputMaximum( (float)calculator->GetMaximum() );
	
	m_PreSigmoidFilter->Update();
	featureImage = m_PreSigmoidFilter->GetOutput();


	
//	m_PreSigmoidFilter->Update();
//	featureImage = m_PreSigmoidFilter->GetOutput();

////////////////////////////////////BEGIN//////////////////////
//	ImageFileWriter<OutputImageType>::Pointer writer = ImageFileWriter<OutputImageType>::New();
	writer->SetInput( featureImage );
	writer->SetFileName("E:\\result\\SPF_feature_After_PreProcess.mha");
	try
	{
	//	writer->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;

	}

///////////////////////////////////////END////////////////////////////


	//获取种子点
	const InputSpatialObjectType * inputSeeds = this->GetInternalInputLandmarks();
	const unsigned int numberOfPoints = inputSeeds->GetNumberOfPoints();
	typedef typename InputSpatialObjectType::PointListType		PointListType;
	const PointListType &points = inputSeeds->GetPoints();

	//存放坐标点
	typedef typename OutputImageType::IndexType					IndexType;
	IndexType  index;

	//存放结果
	outputImage = OutputImageType::New();
	outputImage->SetRegions( featureImage->GetRequestedRegion() );
	outputImage->CopyInformation( featureImage );
	outputImage->Allocate();

	//临时存放生成的梯度值图
	InternalImageType::Pointer  gradientImage = InternalImageType::New();
	gradientImage->SetRegions( featureImage->GetRequestedRegion() );
	gradientImage->CopyInformation( featureImage );
	gradientImage->Allocate();

	//临时存放Signed Pressure Function(SPF)图像
	InternalImageType::Pointer  SPFImage = InternalImageType::New();
	SPFImage->SetRegions( featureImage->GetRequestedRegion() );
	SPFImage->CopyInformation( SPFImage );
	SPFImage->Allocate();

	//初始化featureImage迭代器
	ImageRegionIterator< FeatureImageType > 
		featureIt( featureImage, featureImage->GetRequestedRegion() );

	//初始化SPFImage迭代器
	ImageRegionIterator< InternalImageType > 
		SPFIt( SPFImage, SPFImage->GetRequestedRegion() );

	//初始化outputImage ;背景=-1,初始区域=1.
	ImageRegionIterator< OutputImageType > 
		outputIt( outputImage, outputImage->GetRequestedRegion() );
	for(outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt )
	{
		outputIt.Set(-1.0);
	}
	for( unsigned int i=0; i< numberOfPoints; i++ )
	{
		typename NeighborhoodIterator< OutputImageType >::RadiusType  radius;
		radius.Fill(2);				//初始化区域为到种子点小于7的所有点
		NeighborhoodIterator< OutputImageType >
			outputNeiborIt( radius, outputImage, outputImage->GetRequestedRegion() );

		outputImage->TransformPhysicalPointToIndex( points[i].GetPosition(), index );

		outputNeiborIt.SetLocation( index );

		for(int i=0; i<outputNeiborIt.Size();i++)
		{
			outputNeiborIt.SetPixel(i,1.0);
		}
	}

	////////////////////////////////////BEGIN//////////////////////
	//		ImageFileWriter<OutputImageType>::Pointer writer = ImageFileWriter<OutputImageType>::New();
	writer->SetInput(outputImage);
	writer->SetFileName("E:\\result\\SPF_Initial_Contour.mha");
	try
	{
	//	writer->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;

	}
	///////////////////////////////////////END////////////////////////////

	unsigned int iteratorNum, changedNumber;
	iteratorNum=changedNumber=0;		//初始值=0

	//开始迭代计算，当小于迭代次数并且改变点的数量小于预定值时，继续。
	while( iteratorNum<=m_MaximumNumberOfIterations )
	{	
		//梯度滤波器
		typedef GradientMagnitudeImageFilter< OutputImageType, OutputImageType > GradientFilterType;
		typename GradientFilterType::Pointer m_gradientFilter = GradientFilterType::New();
	//	m_gradientFilter->ReleaseDataFlagOn();

		//高斯滤波器
		typedef DiscreteGaussianImageFilter< InternalImageType, OutputImageType > GaussianFilterType;
		typename GaussianFilterType::Pointer m_gaussianFilter = GaussianFilterType::New();
		m_gaussianFilter->SetUseImageSpacingOn();
		m_gaussianFilter->SetMaximumKernelWidth( m_NeiborSize );
		m_gaussianFilter->SetVariance( m_Sigma );
	//	m_gaussianFilter->ReleaseDataFlagOn();


		m_gradientFilter->SetInput(outputImage);
		m_gradientFilter->Update();
		gradientImage = m_gradientFilter->GetOutput();		//计算梯度幅值图。
////////////////////////////////////BEGIN//////////////////////
		if(iteratorNum%20==0)
		{
	//		ImageFileWriter<OutputImageType>::Pointer writer = ImageFileWriter<OutputImageType>::New();
			writer->SetInput(gradientImage);
			string str="E:\\result\\SPF_gradient_";
			ostringstream oss;
			oss<<iteratorNum;
			str+=oss.str();
			str+=".mha";
			cout<<str<<endl;
			writer->SetFileName(str);
			try
			{
	//			writer->Update();
			}
			catch( itk::ExceptionObject & err )
			{
				std::cerr << "ExceptionObject caught !" << std::endl;
				std::cerr << err << std::endl;

			}
		}
///////////////////////////////////////END////////////////////////////

		//初始化gradientImage迭代器
		ImageRegionIterator< InternalImageType > 
			gradientIt( gradientImage, gradientImage->GetRequestedRegion() );

		float c1,c2;		//用于存放内外均值，c1存放outputImage<0部分,c2存放outputImage>0部分
		c1=c2=0;

		float c1_OnNum, c1_DownNum, c2_OnNum, c2_DownNum;
		c1_OnNum=c1_DownNum=c2_OnNum=c2_DownNum=0;	//分子分母初始化为0

		ImageRegionIterator< OutputImageType > 
			outputIt( outputImage, outputImage->GetRequestedRegion() );

		for( outputIt.GoToBegin(), featureIt.GoToBegin() ; !outputIt.IsAtEnd(); ++outputIt, ++featureIt )
		{
			if( outputIt.Value()<0 )
			{
			//	c1_OnNum+= outputIt.Value()*featureIt.Value();
			//	c1_DownNum+= outputIt.Value();
				c1_OnNum+= featureIt.Value();
				c1_DownNum++;
			}
			else
			{
			//	c2_OnNum+= outputIt.Value()*featureIt.Value();
			//	c2_DownNum+= outputIt.Value();
				c2_OnNum+= featureIt.Value();
				c2_DownNum++;
			}
		}
		c1= c1_OnNum/c1_DownNum;
		c2= c2_OnNum/c2_DownNum;

		float max_SPF = 0;

		for( SPFIt.GoToBegin(), featureIt.GoToBegin() ; !SPFIt.IsAtEnd(); ++SPFIt, ++featureIt )
		{
			SPFIt.Set( featureIt.Value()-(c1+c2)/2 );
			max_SPF = max_SPF>abs(SPFIt.Value())?max_SPF:abs(SPFIt.Value());
		}
		for( gradientIt.GoToBegin(), SPFIt.GoToBegin(), outputIt.GoToBegin(); !gradientIt.IsAtEnd(); ++SPFIt, ++gradientIt, ++outputIt )
		{
			SPFIt.Set( SPFIt.Value()/max_SPF );
	//		cout<<"SPFValue= "<<SPFIt.Value()<<endl;
	//		cout<<"GradientValue="<<gradientIt.Value()<<endl; 
			float temp = m_delta*( m_Propagation*SPFIt.Value()*gradientIt.Value() );
	//		cout<<"Value= "<<outputIt.Value()<<", temp = "<<temp<<endl; 
			outputIt.Set( outputIt.Value() + temp );
	/*		if ( outputIt.Value()<0.0 )
			{
				outputIt.Set( -1.0 );
			}
			else
			{
				outputIt.Set( 1.0 );
			}
			*/
		}

		////////////////////////////////////BEGIN//////////////////////
		
		if(iteratorNum%20==0)
		{
	//		ImageFileWriter<OutputImageType>::Pointer writer = ImageFileWriter<OutputImageType>::New();
			writer->SetInput(outputImage);
			string str="E:\\result\\SPF_OutputImage_Unset_";
			ostringstream oss;
			oss<<iteratorNum;
			str+=oss.str();
			str+=".mha";
			writer->SetFileName(str);
			try
			{
	//			writer->Update();
			}
			catch( itk::ExceptionObject & err )
			{
				std::cerr << "ExceptionObject caught !" << std::endl;
				std::cerr << err << std::endl;

			}
		}
		///////////////////////////////////////END////////////////////////////

/////////////////////////
		for( outputIt.GoToBegin(); !outputIt.IsAtEnd();++outputIt )
		{
			if ( outputIt.Value()<0.0 )
			{
				outputIt.Set( -1.0 );
			}
			else
			{
				outputIt.Set( 1.0 );
			}
		}
//////////////////////////

////////////////////////////////////BEGIN//////////////////////
		if (iteratorNum%20==0)
		{
	//		ImageFileWriter<OutputImageType>::Pointer writer = ImageFileWriter<OutputImageType>::New();
			writer->SetInput(outputImage);
			string str="E:\\result\\SPF_OutputImage_BeforeGaussian_";
			ostringstream oss;
			oss<<iteratorNum;
			str+=oss.str();
			str+=".mha";
			writer->SetFileName(str);
			try
			{
	//			writer->Update();
			}
			catch( itk::ExceptionObject & err )
			{
				std::cerr << "ExceptionObject caught !" << std::endl;
				std::cerr << err << std::endl;

			}
		}
///////////////////////////////////////END////////////////////////////

		m_gaussianFilter->SetInput(outputImage);
		m_gaussianFilter->Update();
		outputImage = m_gaussianFilter->GetOutput();

		iteratorNum++;
////////////////////////////////////BEGIN//////////////////////
		if(iteratorNum%20==0)
		{
	//		ImageFileWriter<OutputImageType>::Pointer writer = ImageFileWriter<OutputImageType>::New();
			writer->SetInput(outputImage);
			string str="E:\\result\\SPF_OutputImage_AfterGaussian_";
			ostringstream oss;
			oss<<iteratorNum;
			str+=oss.str();
			str+=".mha";
			writer->SetFileName("E:\\result\\SPF_OutputImage_AfterGaussian.mha");
			try
			{
	//			writer->Update();
			}
			catch( itk::ExceptionObject & err )
			{
				std::cerr << "ExceptionObject caught !" << std::endl;
				std::cerr << err << std::endl;

			}
		}
///////////////////////////////////////END////////////////////////////

	}

//	m_gradientFilter->SetInput(outputImage);
//	m_gradientFilter->Update();



	
	

////////////////////////////////////BEGIN//////////////////////
//	ImageFileWriter<OutputImageType>::Pointer writer = ImageFileWriter<OutputImageType>::New();
	writer->SetInput(outputImage);
	writer->SetFileName("E:\\result\\SPF_Result_OutputImage.mha");
	try
	{
//		writer->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;

	}
///////////////////////////////////////END////////////////////////////
	this->PackOutputImageInOutputSpatialObject( outputImage );
}

} // end namespace itk

#endif
