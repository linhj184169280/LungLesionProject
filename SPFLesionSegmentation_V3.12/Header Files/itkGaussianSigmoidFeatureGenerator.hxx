/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkGaussianSigmoidFeatureGenerator.hxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGaussianSigmoidFeatureGenerator_hxx
#define __itkGaussianSigmoidFeatureGenerator_hxx

#include "itkGaussianSigmoidFeatureGenerator.h"
#include "itkProgressAccumulator.h"

namespace itk
{

/**
 * Constructor
 */
template <unsigned int NDimension>
GaussianSigmoidFeatureGenerator<NDimension>
::GaussianSigmoidFeatureGenerator()
{
  this->SetNumberOfRequiredInputs( 1 );
  this->SetNumberOfRequiredOutputs( 1 );

  this->m_CastToRealFilter = CastToRealFilterType::New();
  this->m_CastToRealFilter->ReleaseDataFlagOn();

  this->m_SigmoidFilter = SigmoidFilterType::New();
  this->m_SigmoidFilter->ReleaseDataFlagOn();

  this->m_GaussianFilter = GaussianFilterType::New();
  this->m_GaussianFilter->ReleaseDataFlagOn();

  typename OutputImageSpatialObjectType::Pointer outputObject = OutputImageSpatialObjectType::New();

  this->ProcessObject::SetNthOutput( 0, outputObject.GetPointer() );

  this->m_Alpha = -1.0;
  this->m_Beta = 128.0;
  this->m_GaussianSigma = 0.5;
  this->m_GaussianNeiborSize = 5;
}


/*
 * Destructor
 */
template <unsigned int NDimension>
GaussianSigmoidFeatureGenerator<NDimension>
::~GaussianSigmoidFeatureGenerator()
{
}

template <unsigned int NDimension>
void
GaussianSigmoidFeatureGenerator<NDimension>
::SetInput( const SpatialObjectType * spatialObject )
{
  // Process object is not const-correct so the const casting is required.
  this->SetNthInput(0, const_cast<SpatialObjectType *>( spatialObject ));
}

template <unsigned int NDimension>
const typename GaussianSigmoidFeatureGenerator<NDimension>::SpatialObjectType *
GaussianSigmoidFeatureGenerator<NDimension>
::GetFeature() const
{
  return static_cast<const SpatialObjectType*>(this->ProcessObject::GetOutput(0));
}


/*
 * PrintSelf
 */
template <unsigned int NDimension>
void
GaussianSigmoidFeatureGenerator<NDimension>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
}


/*
 * Generate Data
 */
template <unsigned int NDimension>
void
GaussianSigmoidFeatureGenerator<NDimension>
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

  this->m_CastToRealFilter->SetInput( inputImage );

  this->m_GaussianFilter->SetInput( this->m_CastToRealFilter->GetOutput() );
  this->m_GaussianFilter->SetVariance( m_GaussianSigma );
  this->m_GaussianFilter->SetMaximumKernelWidth( m_GaussianNeiborSize );

  this->m_SigmoidFilter->SetInput( this->m_GaussianFilter->GetOutput() );
  this->m_SigmoidFilter->SetAlpha( this->m_Alpha );
  this->m_SigmoidFilter->SetBeta( this->m_Beta );
  this->m_SigmoidFilter->SetOutputMinimum( 0.0 );
  this->m_SigmoidFilter->SetOutputMaximum( 1.0 );

  this->m_SigmoidFilter->Update();

  typename OutputImageType::Pointer outputImage = this->m_SigmoidFilter->GetOutput();

  outputImage->DisconnectPipeline();

  OutputImageSpatialObjectType * outputObject =
    dynamic_cast< OutputImageSpatialObjectType * >(this->ProcessObject::GetOutput(0));

  outputObject->SetImage( outputImage );
}

} // end namespace itk

#endif
