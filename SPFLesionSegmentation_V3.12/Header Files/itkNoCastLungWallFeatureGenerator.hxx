/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkNoCastLungWallFeatureGenerator.hxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkNoCastLungWallFeatureGenerator_hxx
#define __itkNoCastLungWallFeatureGenerator_hxx

#include "itkNoCastLungWallFeatureGenerator.h"
#include "itkProgressAccumulator.h"


namespace itk
{

/**
 * Constructor
 */
template <unsigned int NDimension>
NoCastLungWallFeatureGenerator<NDimension>
::NoCastLungWallFeatureGenerator()
{
  this->SetNumberOfRequiredInputs( 1 );
  this->SetNumberOfRequiredOutputs( 1 );

  this->m_ThresholdFilter = ThresholdFilterType::New();
  this->m_VotingHoleFillingFilter = VotingHoleFillingFilterType::New();

  this->m_ThresholdFilter->ReleaseDataFlagOn();
  this->m_VotingHoleFillingFilter->ReleaseDataFlagOn();

  typename OutputImageSpatialObjectType::Pointer outputObject = OutputImageSpatialObjectType::New();

  this->ProcessObject::SetNthOutput( 0, outputObject.GetPointer() );

  this->m_LungThreshold = -400;
}


/*
 * Destructor
 */
template <unsigned int NDimension>
NoCastLungWallFeatureGenerator<NDimension>
::~NoCastLungWallFeatureGenerator()
{
}

template <unsigned int NDimension>
void
NoCastLungWallFeatureGenerator<NDimension>
::SetInput( const SpatialObjectType * spatialObject )
{
  // Process object is not const-correct so the const casting is required.
  this->SetNthInput(0, const_cast<SpatialObjectType *>( spatialObject ));
}

template <unsigned int NDimension>
const typename NoCastLungWallFeatureGenerator<NDimension>::SpatialObjectType *
NoCastLungWallFeatureGenerator<NDimension>
::GetFeature() const
{
  if (this->GetNumberOfOutputs() < 1)
    {
    return 0;
    }

  return static_cast<const SpatialObjectType*>(this->ProcessObject::GetOutput(0));

}


/*
 * PrintSelf
 */
template <unsigned int NDimension>
void
NoCastLungWallFeatureGenerator<NDimension>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Lung threshold " << this->m_ThresholdFilter << std::endl;
}


/*
 * Generate Data
 */
template <unsigned int NDimension>
void
NoCastLungWallFeatureGenerator<NDimension>
::GenerateData()
{
  typename InputImageSpatialObjectType::ConstPointer inputObject =
    dynamic_cast<const InputImageSpatialObjectType * >( this->ProcessObject::GetInput(0) );

  if( !inputObject )
    {
    itkExceptionMacro("Missing input spatial object");
    }

  InputImageType * inputImage = const_cast< InputImageType * >(inputObject->GetImage());

  if( !inputImage )
    {
    itkExceptionMacro("Missing input image");
    }

  // Report progress.
  ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter(this);
  progress->RegisterInternalFilter( this->m_ThresholdFilter, 0.1 );
  progress->RegisterInternalFilter( this->m_VotingHoleFillingFilter, 0.9 );

  this->m_ThresholdFilter->SetInput( inputImage );
  this->m_VotingHoleFillingFilter->SetInput( this->m_ThresholdFilter->GetOutput() );

  this->m_ThresholdFilter->SetLowerThreshold( this->m_LungThreshold );
  this->m_ThresholdFilter->SetUpperThreshold( 3000 );

  this->m_ThresholdFilter->SetInsideValue( 0.0 );
  this->m_ThresholdFilter->SetOutsideValue( 1.0 );

  typename InternalImageType::SizeType  ballManhattanRadius;

  ballManhattanRadius.Fill( 3 );

  this->m_VotingHoleFillingFilter->SetRadius( ballManhattanRadius );
  this->m_VotingHoleFillingFilter->SetBackgroundValue( 0.0 );
  this->m_VotingHoleFillingFilter->SetForegroundValue( 1.0 );
  this->m_VotingHoleFillingFilter->SetMajorityThreshold( 1 );
  this->m_VotingHoleFillingFilter->SetMaximumNumberOfIterations( 1000 );

  this->m_VotingHoleFillingFilter->Update();

  std::cout << "Used " << this->m_VotingHoleFillingFilter->GetCurrentIterationNumber() << " iterations " << std::endl;
  std::cout << "Changed " << this->m_VotingHoleFillingFilter->GetTotalNumberOfPixelsChanged() << " pixels " << std::endl;

  typename OutputImageType::Pointer internalImage = this->m_VotingHoleFillingFilter->GetOutput();
/////////
//  typename OutputImageType::Pointer outputImage = this->m_VotingHoleFillingFilter->GetOutput();

//  outputImage->DisconnectPipeline();

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
  writer->SetFileName("E:\\result\\LungWall_InternalImage.mha");
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
  writer->SetFileName("E:\\result\\LungWall_Feature_Output.mha");
  try
  {
//	  writer->Update();
  }
  catch( itk::ExceptionObject & err )
  {
	  std::cerr << "ExceptionObject caught !" << std::endl;
	  std::cerr << err << std::endl;

  }

  ImageFileWriter<InputImageType>::Pointer writer2 = ImageFileWriter<InputImageType>::New();
  writer2->SetInput( inputImage );
  writer2->SetFileName("E:\\result\\ROI_LUNGInternalImage.mha");
  try
  {
//	  writer2->Update();
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

} // end namespace itk

#endif
