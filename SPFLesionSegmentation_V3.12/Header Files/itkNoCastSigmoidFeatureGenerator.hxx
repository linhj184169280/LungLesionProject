/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkNoCastSigmoidFeatureGenerator.hxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkNoCastSigmoidFeatureGenerator_hxx
#define __itkNoCastSigmoidFeatureGenerator_hxx

#include "itkNoCastSigmoidFeatureGenerator.h"
#include "itkProgressAccumulator.h"

namespace itk
{

/**
 * Constructor
 */
template <unsigned int NDimension>
NoCastSigmoidFeatureGenerator<NDimension>
::NoCastSigmoidFeatureGenerator()
{
  this->SetNumberOfRequiredInputs( 1 );
  this->SetNumberOfRequiredOutputs( 1 );

  this->m_SigmoidFilter = SigmoidFilterType::New();

  this->m_SigmoidFilter->ReleaseDataFlagOn();

  typename OutputImageSpatialObjectType::Pointer outputObject = OutputImageSpatialObjectType::New();

  this->ProcessObject::SetNthOutput( 0, outputObject.GetPointer() );

  this->m_Alpha = -1.0;
  this->m_Beta = 128.0;
}


/*
 * Destructor
 */
template <unsigned int NDimension>
NoCastSigmoidFeatureGenerator<NDimension>
::~NoCastSigmoidFeatureGenerator()
{
}

template <unsigned int NDimension>
void
NoCastSigmoidFeatureGenerator<NDimension>
::SetInput( const SpatialObjectType * spatialObject )
{
  // Process object is not const-correct so the const casting is required.
  this->SetNthInput(0, const_cast<SpatialObjectType *>( spatialObject ));
}

template <unsigned int NDimension>
const typename NoCastSigmoidFeatureGenerator<NDimension>::SpatialObjectType *
NoCastSigmoidFeatureGenerator<NDimension>
::GetFeature() const
{
  return static_cast<const SpatialObjectType*>(this->ProcessObject::GetOutput(0));
}


/*
 * PrintSelf
 */
template <unsigned int NDimension>
void
NoCastSigmoidFeatureGenerator<NDimension>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
}


/*
 * Generate Data
 */
template <unsigned int NDimension>
void
NoCastSigmoidFeatureGenerator<NDimension>
::GenerateData()
{
  // Report progress.
  ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter(this);
  progress->RegisterInternalFilter( this->m_SigmoidFilter, 1.0 );  

  typename InputImageSpatialObjectType::ConstPointer inputObject =
    dynamic_cast<const InputImageSpatialObjectType * >( this->ProcessObject::GetInput(0) );

  if( !inputObject )
    {
    itkExceptionMacro("Missing input spatial object or incorrect type");
    }

  const InputImageType * inputImage = inputObject->GetImage();

  if( !inputImage )
    {
    itkExceptionMacro("Missing input image");
    }
/////////////////////////////////////////////////////////
  ImageFileWriter<InputImageType>::Pointer writer2 = ImageFileWriter<InputImageType>::New();
  writer2->SetInput( inputImage );
  writer2->SetFileName("E:\\result\\ROI_SigmoidImage.mha");
  try
  {
//	  writer2->Update();
  }
  catch( itk::ExceptionObject & err )
  {
	  std::cerr << "ExceptionObject caught !" << std::endl;
	  std::cerr << err << std::endl;

  }
///////////////////////////////////////////////////

  typedef MinimumMaximumImageCalculator< InputImageType > CalculatorType;
  typename CalculatorType::Pointer calculator = CalculatorType::New();
  calculator->SetImage( inputImage );
  calculator->Compute();

  this->m_SigmoidFilter->SetInput( inputImage );
  this->m_SigmoidFilter->SetAlpha( this->m_Alpha );
  this->m_SigmoidFilter->SetBeta( this->m_Beta );
  this->m_SigmoidFilter->SetOutputMinimum( (float)calculator->GetMinimum() );
  this->m_SigmoidFilter->SetOutputMaximum( (float)calculator->GetMaximum() );
//  this->m_SigmoidFilter->SetOutputMinimum( 0.0 );
//  this->m_SigmoidFilter->SetOutputMaximum( 1.0 );

  cout<<"Sigmoid: Alpha="<<m_Alpha<<"  Beta="<<m_Beta<<endl;

  this->m_SigmoidFilter->Update();

  typename OutputImageType::Pointer outputImage = this->m_SigmoidFilter->GetOutput();

  outputImage->DisconnectPipeline();

  OutputImageSpatialObjectType * outputObject =
    dynamic_cast< OutputImageSpatialObjectType * >(this->ProcessObject::GetOutput(0));

  outputObject->SetImage( outputImage );
}

} // end namespace itk

#endif
