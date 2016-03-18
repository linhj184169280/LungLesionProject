/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkNoCastCannyEdgesFeatureGenerator.hxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkNoCastCannyEdgesFeatureGenerator_hxx
#define __itkNoCastCannyEdgesFeatureGenerator_hxx

#include "itkNoCastCannyEdgesFeatureGenerator.h"


namespace itk
{

/**
 * Constructor
 */
template <unsigned int NDimension>
NoCastCannyEdgesFeatureGenerator<NDimension>
::NoCastCannyEdgesFeatureGenerator()
{
  this->SetNumberOfRequiredInputs( 1 );

  this->m_CastFilter        = CastFilterType::New();			//用于把int型图像转换成float型
  this->m_RescaleFilter     = RescaleFilterType::New();
  this->m_CannyFilter       = CannyEdgeFilterType::New();

  typename OutputImageSpatialObjectType::Pointer 
    outputObject = OutputImageSpatialObjectType::New();

  this->ProcessObject::SetNthOutput( 0, outputObject.GetPointer() );

  this->m_Sigma.Fill(1.0);
  this->m_UpperThreshold = NumericTraits< InternalPixelType >::max();
  this->m_LowerThreshold = NumericTraits< InternalPixelType >::min();

  this->m_RescaleFilter->SetOutputMinimum( 1.0 );
  this->m_RescaleFilter->SetOutputMaximum( 0.0 );

  this->m_RescaleFilter->SetWindowMinimum( 0.0 );
  this->m_RescaleFilter->SetWindowMaximum( 1.0 );
}


/*
 * Destructor
 */
template <unsigned int NDimension>
NoCastCannyEdgesFeatureGenerator<NDimension>
::~NoCastCannyEdgesFeatureGenerator()
{
}

template <unsigned int NDimension>
void
NoCastCannyEdgesFeatureGenerator<NDimension>
::SetInput( const SpatialObjectType * spatialObject )
{
  // Process object is not const-correct so the const casting is required.
  this->SetNthInput(0, const_cast<SpatialObjectType *>( spatialObject ));
}

template <unsigned int NDimension>
const typename NoCastCannyEdgesFeatureGenerator<NDimension>::SpatialObjectType *
NoCastCannyEdgesFeatureGenerator<NDimension>
::GetFeature() const
{
  return static_cast<const SpatialObjectType*>(this->ProcessObject::GetOutput(0));
}


/*
 * PrintSelf
 */
template <unsigned int NDimension>
void
NoCastCannyEdgesFeatureGenerator<NDimension>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
}


/*
 * Generate Data
 */
template <unsigned int NDimension>
void
NoCastCannyEdgesFeatureGenerator<NDimension>
::GenerateData()
{
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

  this->m_CastFilter->SetInput( inputImage );
  this->m_CannyFilter->SetInput( this->m_CastFilter->GetOutput() );
  this->m_RescaleFilter->SetInput( this->m_CannyFilter->GetOutput() );

  this->m_CannyFilter->SetSigmaArray( this->m_Sigma );
  this->m_CannyFilter->SetUpperThreshold( this->m_UpperThreshold );
  this->m_CannyFilter->SetLowerThreshold( this->m_LowerThreshold );
  this->m_CannyFilter->SetOutsideValue(NumericTraits<InternalPixelType>::Zero);

  this->m_RescaleFilter->Update();
//////
 // typename OutputImageType::Pointer outputImage = this->m_RescaleFilter->GetOutput();
//  outputImage->DisconnectPipeline();

  typename OutputImageType::Pointer internalImage = this->m_RescaleFilter->GetOutput();

  internalImage->DisconnectPipeline();

  typename OutputImageType::Pointer outputImage = OutputImageType::New();
  outputImage->SetRegions( internalImage->GetRequestedRegion() );
  outputImage->CopyInformation( internalImage );
  outputImage->Allocate();

  ImageRegionIterator< InputImageType > 
	  inputIt( inputImage, inputImage->GetRequestedRegion() );

  ImageRegionIterator< OutputImageType > 
	  internalIt( internalImage, internalImage->GetRequestedRegion() );

  ImageRegionIterator< OutputImageType > 
	  outputIt( outputImage, outputImage->GetRequestedRegion() );

  typedef MinimumMaximumImageCalculator< InputImageType > CalculatorType;
  typename CalculatorType::Pointer calculator = CalculatorType::New();
  calculator->SetImage( inputImage );
  calculator->Compute();
  float minValue = (float)calculator->GetMinimum();
  float maxValue = (float)calculator->GetMaximum();
//  float minValue = 0.0;
//  float maxValue = 1.0;

  for (inputIt.GoToBegin(),internalIt.GoToBegin(),outputIt.GoToBegin();!internalIt.IsAtEnd();++inputIt,++internalIt,++outputIt)
  {
	  if (internalIt.Value()<0.5)
	  {
		  outputIt.Set( minValue );
	  }
	  else
	  {
		  outputIt.Set( (float)inputIt.Value() );
		//  outputIt.Set( maxValue );
	  }

  }

  ////////////////////////////////////BEGIN//////////////////////
  ImageFileWriter<OutputImageType>::Pointer writer = ImageFileWriter<OutputImageType>::New();
  writer->SetInput( internalImage );
  writer->SetFileName("E:\\result\\Canny_InternalImage.mha");
  try
  {
//	  writer->Update();
  }
  catch( itk::ExceptionObject & err )
  {
	  std::cerr << "ExceptionObject caught !" << std::endl;
	  std::cerr << err << std::endl;

  }

  writer->SetInput( outputImage );
  writer->SetFileName("E:\\result\\Canny_Feature_Output.mha");
  try
  {
//	  writer->Update();
  }
  catch( itk::ExceptionObject & err )
  {
	  std::cerr << "ExceptionObject caught !" << std::endl;
	  std::cerr << err << std::endl;

  }

  ///////////////////////////////////////END////////////////////////////
  

  OutputImageSpatialObjectType * outputObject =
    dynamic_cast< OutputImageSpatialObjectType * >(this->ProcessObject::GetOutput(0));

  outputObject->SetImage( outputImage );
}


// Set value of Sigma (isotropic)

template <unsigned int NDimension>
void 
NoCastCannyEdgesFeatureGenerator<NDimension>
::SetSigma( ScalarRealType sigma )
{
  SigmaArrayType sigmas(sigma);
  this->SetSigmaArray(sigmas);
}


// Set value of Sigma (an-isotropic)

template <unsigned int NDimension>
void 
NoCastCannyEdgesFeatureGenerator<NDimension>
::SetSigmaArray( const SigmaArrayType & sigma )
{
  if (this->m_Sigma != sigma)
    {
    this->m_Sigma = sigma;
    this->Modified();
    }
}


// Get the sigma array.
template <unsigned int NDimension>
typename NoCastCannyEdgesFeatureGenerator<NDimension>::SigmaArrayType
NoCastCannyEdgesFeatureGenerator<NDimension>
::GetSigmaArray() const
{
  return m_Sigma;
}


// Get the sigma scalar. If the sigma is anisotropic, we will just
// return the sigma along the first dimension.
template <unsigned int NDimension>
typename NoCastCannyEdgesFeatureGenerator<NDimension>::ScalarRealType
NoCastCannyEdgesFeatureGenerator<NDimension>
::GetSigma() const
{
  return m_Sigma[0];
}

} // end namespace itk

#endif
