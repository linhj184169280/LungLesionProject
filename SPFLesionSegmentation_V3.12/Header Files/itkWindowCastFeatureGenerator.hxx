/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkWindowCastFeatureGenerator.hxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkWindowCastFeatureGenerator_hxx
#define __itkWindowCastFeatureGenerator_hxx

#include "itkWindowCastFeatureGenerator.h"
#include "itkProgressAccumulator.h"

namespace itk
{

/**
 * Constructor
 */
template <unsigned int NDimension>
WindowCastFeatureGenerator<NDimension>
::WindowCastFeatureGenerator()
{
  this->SetNumberOfRequiredInputs( 1 );
  this->SetNumberOfRequiredOutputs( 1 );

  this->m_WindowCastFilter = WindowCastFilterType::New();

  this->m_WindowCastFilter->ReleaseDataFlagOn();

  typename OutputImageSpatialObjectType::Pointer outputObject = OutputImageSpatialObjectType::New();

  this->ProcessObject::SetNthOutput( 0, outputObject.GetPointer() );

  this->m_OutputMinimum = 0.0;
  this->m_OutputMaximum = 1.0;
}


/*
 * Destructor
 */
template <unsigned int NDimension>
WindowCastFeatureGenerator<NDimension>
::~WindowCastFeatureGenerator()
{
}

template <unsigned int NDimension>
void
WindowCastFeatureGenerator<NDimension>
::SetInput( const SpatialObjectType * spatialObject )
{
  // Process object is not const-correct so the const casting is required.
  this->SetNthInput(0, const_cast<SpatialObjectType *>( spatialObject ));
}

template <unsigned int NDimension>
const typename WindowCastFeatureGenerator<NDimension>::SpatialObjectType *
WindowCastFeatureGenerator<NDimension>
::GetFeature() const
{
  return static_cast<const SpatialObjectType*>(this->ProcessObject::GetOutput(0));
}


/*
 * PrintSelf
 */
template <unsigned int NDimension>
void
WindowCastFeatureGenerator<NDimension>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
}


/*
 * Generate Data
 */
template <unsigned int NDimension>
void
WindowCastFeatureGenerator<NDimension>
::GenerateData()
{
  // Report progress.
  ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter(this);
  progress->RegisterInternalFilter( this->m_WindowCastFilter, 1.0 );  

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

  typename CalculatorType::Pointer calculator = CalculatorType::New();
  calculator->SetImage( inputImage );
  calculator->Compute();

  this->m_WindowCastFilter->SetInput( inputImage );
  this->m_WindowCastFilter->SetWindowMinimum( calculator->GetMinimum() );
  this->m_WindowCastFilter->SetWindowMaximum( calculator->GetMaximum() ); 
  this->m_WindowCastFilter->SetOutputMinimum( 0.0 );
  this->m_WindowCastFilter->SetOutputMaximum( 1.0 );
  this->m_WindowCastFilter->InPlaceOn();
  this->m_WindowCastFilter->Update();

  typename OutputImageType::Pointer outputImage = this->m_WindowCastFilter->GetOutput();

  outputImage->DisconnectPipeline();

  OutputImageSpatialObjectType * outputObject =
    dynamic_cast< OutputImageSpatialObjectType * >(this->ProcessObject::GetOutput(0));

  outputObject->SetImage( outputImage );
}

} // end namespace itk

#endif
