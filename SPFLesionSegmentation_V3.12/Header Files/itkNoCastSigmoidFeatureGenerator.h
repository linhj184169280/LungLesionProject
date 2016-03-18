/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkNoCastSigmoidFeatureGenerator.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkNoCastSigmoidFeatureGenerator_h
#define __itkNoCastSigmoidFeatureGenerator_h

#include "itkFeatureGenerator.h"
#include "itkImage.h"
#include "itkImageSpatialObject.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkSigmoidImageFilter.h"

namespace itk
{

/** \class NoCastSigmoidFeatureGenerator
 * \brief Generates a feature image by computing the Sigmoid intensity transformation.
 *
 * The typical use of this class would be to generate the feature image needed
 * by a Level Set filter to internally compute its speed image. This
 * transformation is very close to a simply thresholding selection on the input
 * image, but with the advantage of a smooth transition of intensities.
 *
 * SpatialObjects are used as inputs and outputs of this class.
 *
 * \ingroup SpatialObjectFilters
 * \ingroup ITKLesionSizingToolkit
 */
template <unsigned int NDimension>
class ITK_EXPORT NoCastSigmoidFeatureGenerator : public FeatureGenerator<NDimension>
{
public:
  /** Standard class typedefs. */
  typedef NoCastSigmoidFeatureGenerator          Self;
  typedef FeatureGenerator<NDimension>     Superclass;
  typedef SmartPointer<Self>               Pointer;
  typedef SmartPointer<const Self>         ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(NoCastSigmoidFeatureGenerator, FeatureGenerator);

  /** Dimension of the space */
  itkStaticConstMacro(Dimension, unsigned int, NDimension);

  /** Type of spatialObject that will be passed as input to this
   * feature generator. */
  typedef signed short                                      InputPixelType;
  typedef Image< InputPixelType, Dimension >                InputImageType;
  typedef ImageSpatialObject< NDimension, InputPixelType >  InputImageSpatialObjectType;
  typedef typename InputImageSpatialObjectType::Pointer     InputImageSpatialObjectPointer;
  typedef typename Superclass::SpatialObjectType            SpatialObjectType;

  /** Input data that will be used for generating the feature. */
  using ProcessObject::SetInput;
  void SetInput( const SpatialObjectType * input );
  const SpatialObjectType * GetInput() const;

  /** Output data that carries the feature in the form of a
   * SpatialObject. */
  const SpatialObjectType * GetFeature() const;

  /** Alpha value to be used in the Sigmoid filter. */
  itkSetMacro( Alpha, double );
  itkGetMacro( Alpha, double );

  /** Beta value to be used in the Sigmoid filter. */
  itkSetMacro( Beta, double );
  itkGetMacro( Beta, double );

protected:
  NoCastSigmoidFeatureGenerator();
  virtual ~NoCastSigmoidFeatureGenerator();
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Method invoked by the pipeline in order to trigger the computation of
   * the segmentation. */
  void  GenerateData ();

private:
  NoCastSigmoidFeatureGenerator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  typedef float                                       OutputPixelType;
  typedef Image< OutputPixelType, Dimension >         OutputImageType;

  typedef ImageSpatialObject< NDimension, OutputPixelType >  OutputImageSpatialObjectType;

  typedef SigmoidImageFilter<
    InputImageType, OutputImageType >                 SigmoidFilterType;
  typedef typename SigmoidFilterType::Pointer         SigmoidFilterPointer;



  SigmoidFilterPointer            m_SigmoidFilter;

  double                          m_Alpha;
  double                          m_Beta;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
# include "itkNoCastSigmoidFeatureGenerator.hxx"
#endif

#endif
