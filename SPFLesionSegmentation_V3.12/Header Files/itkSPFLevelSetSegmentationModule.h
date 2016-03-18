/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkSPFLevelSetSegmentationModule.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkSPFLevelSetSegmentationModule_h
#define __itkSPFLevelSetSegmentationModule_h

#include "itkSegmentationModule.h"
#include "itkImageSpatialObject.h"
#include "itkLandmarkSpatialObject.h"

namespace itk
{

/** \class SPFLevelSetSegmentationModule
 * \brief Class applies a single-phase level set segmentation method
 *
 * SpatialObjects are used as inputs and outputs of this class.
 *
 * \ingroup SpatialObjectFilters
 * \ingroup ITKLesionSizingToolkit
 */
template <unsigned int NDimension>
class ITK_EXPORT SPFLevelSetSegmentationModule : public SegmentationModule<NDimension>
{
public:
  /** Standard class typedefs. */
  typedef SPFLevelSetSegmentationModule				    Self;
  typedef SegmentationModule<NDimension>                Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SPFLevelSetSegmentationModule, SegmentationModule);

  /** Dimension of the space */
  itkStaticConstMacro(Dimension, unsigned int, NDimension);

  /** Type of spatialObject that will be passed as input and output of this
   * segmentation method. */
  typedef typename Superclass::SpatialObjectType         SpatialObjectType;
  typedef typename Superclass::SpatialObjectPointer      SpatialObjectPointer;

  /** Types of the input, feature and output images. */
  typedef float                                         InputImageType;
  typedef float                                         InternalPixelType;
  typedef float                                         FeaturePixelType;
  typedef float                                         OutputPixelType;
  typedef Image< InternalPixelType, NDimension >        InternalImageType;
  typedef Image< FeaturePixelType, NDimension >         FeatureImageType;
  typedef Image< OutputPixelType, NDimension >          OutputImageType;

  /** Types of the Spatial objects used for input, feature and output images. */
  typedef LandmarkSpatialObject< NDimension >     InputSpatialObjectType;
  typedef ImageSpatialObject< NDimension, FeaturePixelType >   FeatureSpatialObjectType;
  typedef ImageSpatialObject< NDimension, OutputPixelType >    OutputSpatialObjectType;

  

  /** Gaussian sigma. */
  itkSetMacro( Sigma, double );
  itkGetMacro( Sigma, double );

  /** Gaussian NeiborSize */
  itkSetMacro( NeiborSize, unsigned int );
  itkGetMacro( NeiborSize, unsigned int );

  /** Parameter that controls the updata : u = u + m_delta*(increase value) . */
  itkSetMacro( delta, double );
  itkGetMacro( delta, double );

  /** Gaussian sigma. */
  itkSetMacro( Propagation, double );
  itkGetMacro( Propagation, double );

  /** when the number of the points changed lower than the MinChangeNumber ,then the Iteration will stop. */
  itkSetMacro( MinChangeNumber, unsigned int );
  itkGetMacro( MinChangeNumber, unsigned int );

  /** Maximum number of iterations that the level set solve will run. */
  itkSetMacro( MaximumNumberOfIterations, unsigned int );
  itkGetMacro( MaximumNumberOfIterations, unsigned int );

  /** Invert the output image. This is a convenience method intended to make
   * uniform the convention that segmentations are encoded with positive values
   * in the pixels inside of the segmented object and negative values in the
   * pixels outside of the segmented object. This is opposed to the general
   * convention of ITK level sets, where the values inside the object are
   * negative, and for this reason they must be inverted here. By default the
   * intensities must be inverted, and therefore, by default this variable will
   * be set to true. However, when combining multiple level sets in a sequence,
   * this variable should be set to false. */
  itkSetMacro( InvertOutputIntensities, bool );
  itkGetMacro( InvertOutputIntensities, bool );
  itkBooleanMacro( InvertOutputIntensities );
  
protected:
  SPFLevelSetSegmentationModule();
  virtual ~SPFLevelSetSegmentationModule();
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Method invoked by the pipeline in order to trigger the computation of
   * the segmentation. */
  void  GenerateData ();

  /** Set the output image as cargo of the output SpatialObject. */
  void PackOutputImageInOutputSpatialObject( OutputImageType * outputImage );

  /** Extract the input set of landmark points to be used as seeds. */
  const InputSpatialObjectType * GetInternalInputLandmarks() const;

  /** Extract the input feature image from the input feature spatial object. */
  const FeatureImageType * GetInternalFeatureImage() const;

  /** Using to store the output */
  typename OutputImageType::Pointer outputImage;

private:
  SPFLevelSetSegmentationModule(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  double			m_Sigma;
  unsigned short    m_NeiborSize;
  double			m_delta;
  double			m_Propagation;

  unsigned int  m_MaximumNumberOfIterations;
  unsigned int  m_MinChangeNumber;

  bool          m_InvertOutputIntensities;

//  typedef typename InputImageType::ConstPointer  ImageConstPointer;
//  mutable ImageConstPointer m_ZeroSetInputImage;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
# include "itkSPFLevelSetSegmentationModule.hxx"
#endif

#endif
