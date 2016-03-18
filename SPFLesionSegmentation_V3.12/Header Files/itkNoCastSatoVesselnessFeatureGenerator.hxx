/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkNoCastSatoVesselnessFeatureGenerator.hxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkNoCastSatoVesselnessFeatureGenerator_hxx
#define __itkNoCastSatoVesselnessFeatureGenerator_hxx

#include "itkNoCastSatoVesselnessFeatureGenerator.h"
#include "itkProgressAccumulator.h"


namespace itk
{

/**
 * Constructor
 */
template <unsigned int NDimension>
NoCastSatoVesselnessFeatureGenerator<NDimension>
::NoCastSatoVesselnessFeatureGenerator()
{
  this->m_SigmoidFilter = SigmoidFilterType::New();

  this->m_SigmoidFilter->ReleaseDataFlagOn();

  this->m_SigmoidAlpha =  -1.0;
  this->m_SigmoidBeta = 90.0;
}


/*
 * Destructor
 */
template <unsigned int NDimension>
NoCastSatoVesselnessFeatureGenerator<NDimension>
::~NoCastSatoVesselnessFeatureGenerator()
{
}


/**
 * PrintSelf
 */
template <unsigned int NDimension>
void
NoCastSatoVesselnessFeatureGenerator<NDimension>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Sigmoid Alpha " << this->m_SigmoidAlpha << std::endl;
  os << indent << "Sigmoid Beta " << this->m_SigmoidBeta << std::endl;
}


/*
 * Generate Data
 */
template <unsigned int NDimension>
void
NoCastSatoVesselnessFeatureGenerator<NDimension>
::GenerateData()
{
	this->Superclass::GenerateData();

	// Report progress. Actually, the superclass will report upto 1 in
	// the superclass's generate data method. This will start again
	// from 0, but that's ok. :)
	ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
	progress->SetMiniPipelineFilter(this);
	progress->RegisterInternalFilter( this->m_SigmoidFilter, 1.0 ); 

	//
	// Take the output of the superclass, and do further processing on it.
	// 
	typename OutputImageSpatialObjectType::Pointer outputObject =
		dynamic_cast<OutputImageSpatialObjectType * >( this->GetInternalFeature() );

	OutputImageType * internalImage = const_cast< OutputImageType * >( outputObject->GetImage() );

	////获取未经处理的源图像
	typename InputImageSpatialObjectType::ConstPointer inputObject =
		dynamic_cast<const InputImageSpatialObjectType * >( this->ProcessObject::GetInput(0) );

	if( !inputObject )
	{
		itkExceptionMacro("Missing input spatial object or incorrect type");
	}

	InputImageType * inputImage = const_cast< InputImageType * >( inputObject->GetImage() );

	if( !inputImage )
	{
		itkExceptionMacro("Missing input image");
	}

	this->m_SigmoidFilter->SetInput( internalImage );

	this->m_SigmoidFilter->SetAlpha( this->m_SigmoidAlpha );
	this->m_SigmoidFilter->SetBeta( this->m_SigmoidBeta );

	typedef MinimumMaximumImageCalculator< InputImageType > CalculatorType;
	typename CalculatorType::Pointer calculator = CalculatorType::New();
	calculator->SetImage( inputImage );
	calculator->Compute();

	float minValue = (float)calculator->GetMinimum();
	float maxValue = (float)calculator->GetMaximum();

	this->m_SigmoidFilter->SetOutputMinimum( minValue );
	this->m_SigmoidFilter->SetOutputMaximum( maxValue );
//	this->m_SigmoidFilter->SetOutputMinimum( 0.0 );
//	this->m_SigmoidFilter->SetOutputMaximum( 1.0 );

	this->m_SigmoidFilter->Update();

//	internalImage = this->m_SigmoidFilter->GetOutput();

	typename OutputImageType::Pointer outputImage = this->m_SigmoidFilter->GetOutput();
	outputImage->DisconnectPipeline();
//	outputObject->SetImage( outputImage );
/*
	typename OutputImageType::Pointer outputImage = OutputImageType::New();
	outputImage->SetRegions( internalImage->GetRequestedRegion() );
	outputImage->CopyInformation( internalImage );
	outputImage->Allocate();

//	outputImage = internalImage;
	ImageRegionIterator< InputImageType > 
		inputIt( inputImage, inputImage->GetRequestedRegion() );

	ImageRegionIterator< OutputImageType > 
		internalIt( internalImage, internalImage->GetRequestedRegion() );

	ImageRegionIterator< OutputImageType > 
		outputIt( outputImage, outputImage->GetRequestedRegion() );

	for (inputIt.GoToBegin(),internalIt.GoToBegin(),outputIt.GoToBegin();!inputIt.IsAtEnd();++inputIt,++internalIt,++outputIt)
	{
		if (internalIt.Value()<0.1)
		{
			outputIt.Set( minValue );
		}
		else
		{
			outputIt.Set( (float)inputIt.Value() );
		}
		//outputIt.Set( internalIt.Value() );
	}

	////////////////////////////////////BEGIN//////////////////////
	ImageFileWriter<OutputImageType>::Pointer writer = ImageFileWriter<OutputImageType>::New();
	writer->SetInput( internalImage );
	writer->SetFileName("E:\\result\\Vessel_InternalImage.mha");
	try
	{
		writer->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;

	}

	writer->SetInput( outputImage );
	writer->SetFileName("E:\\result\\Vessel_Feature_Output.mha");
	try
	{
		writer->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;

	}
*/
	///////////////////////////////////////END////////////////////////////

//	outputImage->DisconnectPipeline();

//	internalImage->DisconnectPipeline();
	outputObject->SetImage( outputImage );
	
}

} // end namespace itk

#endif
