/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    FMandGACLevelSetSegmentationModule.hxx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkFMandGACLevelSetSegmentationModule_hxx
#define __itkFMandGACLevelSetSegmentationModule_hxx

#include "itkFMandGACLevelSetSegmentationModule.h"
#include "itkGeodesicActiveContourLevelSetImageFilter.h"
#include "itkProgressAccumulator.h"
#include "itkSigmoidImageFilter.h"
#include "itkNeighborhoodAlgorithm.h"

namespace itk
{


/**
 * Constructor
 */
template <unsigned int NDimension>
FMandGACLevelSetSegmentationModule<NDimension>
::FMandGACLevelSetSegmentationModule()
{
  this->m_FastMarchingModule = FastMarchingModuleType::New();
  this->m_FastMarchingModule->SetDistanceFromSeeds(1.0);
  this->m_FastMarchingModule->SetStoppingValue( 100.0 );
  this->m_FastMarchingModule->InvertOutputIntensitiesOff();
  this->m_GeodesicActiveContourLevelSetModule = GeodesicActiveContourLevelSetModuleType::New();
  this->m_GeodesicActiveContourLevelSetModule->InvertOutputIntensitiesOff();
  this->m_SPFLevelSetModule = SPFLevelSetModuleType::New();
}


/**
 * Destructor
 */
template <unsigned int NDimension>
FMandGACLevelSetSegmentationModule<NDimension>
::~FMandGACLevelSetSegmentationModule()
{
}


/**
 * PrintSelf
 */
template <unsigned int NDimension>
void
FMandGACLevelSetSegmentationModule<NDimension>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
}


/**
 * Generate Data
 */
template <unsigned int NDimension>
void
FMandGACLevelSetSegmentationModule<NDimension>
::GenerateData()
{
  this->m_SPFLevelSetModule->SetInput( this->GetInput() );
  this->m_SPFLevelSetModule->SetFeature( this->GetFeature() );
  this->m_SPFLevelSetModule->Update();

  typedef SigmoidImageFilter<
	  InputImageType, OutputImageType >                 SigmoidFilterType;
  typedef typename SigmoidFilterType::Pointer         SigmoidFilterPointer;
  SigmoidFilterPointer            m_SigmoidFilter;
  m_SigmoidFilter = SigmoidFilterType::New();
//  m_SigmoidFilter->ReleaseDataFlagOn();

  FeatureSpatialObjectType * featureObject =
	  const_cast< FeatureSpatialObjectType * >(dynamic_cast< const FeatureSpatialObjectType * >( this->GetFeature() ));
  FeatureImageType * featureImage = const_cast< FeatureImageType *  >(featureObject->GetImage());

  m_SigmoidFilter->SetInput( featureImage );
  m_SigmoidFilter->SetAlpha( 100 );
  m_SigmoidFilter->SetBeta( -200 );
  m_SigmoidFilter->SetOutputMinimum( 0.0 );
  m_SigmoidFilter->SetOutputMaximum( 1.0 );
  m_SigmoidFilter->Update();
  
  OutputSpatialObjectType::Pointer outputObject = OutputSpatialObjectType::New();
  outputObject->CopyInformation( m_SPFLevelSetModule->GetOutput() );
  outputObject->SetImage( dynamic_cast< const OutputSpatialObjectType * >(m_SPFLevelSetModule->GetOutput())->GetImage() );

  FeatureSpatialObjectType::Pointer featureObject2 = FeatureSpatialObjectType::New();
  featureObject2->CopyInformation( this->GetFeature() );
  featureObject2->SetImage( m_SigmoidFilter->GetOutput() );

  m_GeodesicActiveContourLevelSetModule->SetInput( outputObject );
  m_GeodesicActiveContourLevelSetModule->SetFeature( featureObject2 );
  m_GeodesicActiveContourLevelSetModule->SetMaximumRMSError( this->GetMaximumRMSError() );
  m_GeodesicActiveContourLevelSetModule->SetMaximumNumberOfIterations( this->GetMaximumNumberOfIterations() );
  m_GeodesicActiveContourLevelSetModule->SetPropagationScaling( this->GetPropagationScaling() );
  m_GeodesicActiveContourLevelSetModule->SetCurvatureScaling( this->GetCurvatureScaling() );
  m_GeodesicActiveContourLevelSetModule->SetAdvectionScaling( this->GetAdvectionScaling() );
  m_GeodesicActiveContourLevelSetModule->Update();
////////////////////////
 OutputImageType * outputImage = const_cast< OutputImageType * >(
	  dynamic_cast< const OutputSpatialObjectType * >(
	  m_GeodesicActiveContourLevelSetModule->GetOutput())->GetImage());

/*   ImageRegionIterator< OutputImageType > 
	  outputIt( outputImage, outputImage->GetRequestedRegion() );

  for(outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt )
  {
	  if ( outputIt.Value()<=0.5 )
	  {
		  outputIt.Set( -4.0 );
	  }
	  else
	  {
		  outputIt.Set( 4.0 );
	  }
  }

  this->PackOutputImageInOutputSpatialObject( outputImage );
 */ 
///////////////////////
  ////////////////////////////////////BEGIN//////////////////////
  ImageFileWriter<OutputImageType>::Pointer writer = ImageFileWriter<OutputImageType>::New();
  writer->SetInput( outputImage );
  writer->SetFileName("E:\\result\\SPFGACOutput.mha");
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

  
  this->PackOutputImageInOutputSpatialObject( const_cast< OutputImageType * >(
        dynamic_cast< const OutputSpatialObjectType * >(
        m_GeodesicActiveContourLevelSetModule->GetOutput())->GetImage()) );

}

} // end namespace itk

#endif
