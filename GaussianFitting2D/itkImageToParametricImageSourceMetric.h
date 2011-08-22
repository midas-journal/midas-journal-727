/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkImageToParametricImageSourceMetric.h,v $
  Language:  C++
  Date:      $Date: 2010/04/19 18:50:02 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkImageToParametricImageSourceMetric_h
#define __itkImageToParametricImageSourceMetric_h

// First make sure that the configuration is available.
// This line can be removed once the optimized versions
// gets integrated into the main directories.
#include "itkConfigure.h"

#include "itkImageToImageMetric.h"

namespace itk
{
  
/** \class ImageToParametricImageSourceMetric
 * \brief Computes similarity between two images, one of which is fixed and
 * the other generated from a moving ParametricImageSource.
 *
 * This class computes a value that measures the similarity 
 * between the Fixed image and the parametric Moving image. "Moving"
 * in this metric and subclasses refers to changes in the parameters used
 * to generate the image, not movement induced by a spatial transformation.
 *
 * To compute the similarity, this class uses a delegate ImageToImageMetric
 * that needs to be set with SetImageToImageMetric method prior to using
 * this class.
 *
 * This class is parameterized over two types. The first template argument
 * is the fixed image data and the second template argument is the source of
 * the moving ParametricImageSource.
 *
 * This class was developed by Cory Quammen at the Center for Computer
 * Integrated Systems for Microscopy and Manipulation (http:://www.cismm.org)
 * at the University of North Carolina at Chapel Hill and was supported by
 * NIH grant P41-EB002025-25A1.
 *
 * \ingroup RegistrationMetrics
 *
 */

template <class TFixedImage,  class TMovingImageSource> 
class ITK_EXPORT ImageToParametricImageSourceMetric : public SingleValuedCostFunction 
{
public:
  /** Standard class typedefs. */
  typedef ImageToParametricImageSourceMetric Self;
  typedef SingleValuedCostFunction           Superclass;
  typedef SmartPointer<Self>                 Pointer;
  typedef SmartPointer<const Self>           ConstPointer;

  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageToParametricImageSourceMetric, SingleValuedCostFunction);

  /**  Type of the moving Image. */
  typedef TMovingImageSource 
    MovingImageSourceType;
  typedef typename MovingImageSourceType::OutputImagePixelType
    MovingImageSourcePixelType;
  typedef typename MovingImageSourceType::Pointer
    MovingImageSourcePointer;
  typedef typename MovingImageSourceType::OutputImageType
    MovingImageSourceOutputImageType;
  typedef typename MovingImageSourceOutputImageType::Pointer
    MovingImageSourceOutputImagePointerType;

  /**  Type of the fixed Image. */
  typedef TFixedImage                                FixedImageType;
  typedef typename FixedImageType::ConstPointer      FixedImageConstPointer;
  typedef typename FixedImageType::RegionType        FixedImageRegionType;

  /**  Type of the delegate image comparison metric. */
  typedef ImageToImageMetric<FixedImageType, MovingImageSourceOutputImageType>
    DelegateMetricType;
  typedef typename DelegateMetricType::Pointer
    DelegateMetricTypePointer;

  /** Constants for the image dimensions */
  itkStaticConstMacro(MovingImageSourceDimension, 
                      unsigned int,
                      MovingImageSourceType::OutputImageDimension);
  itkStaticConstMacro(FixedImageDimension, 
                      unsigned int,
                      FixedImageType::ImageDimension);
  
  typedef typename NumericTraits<MovingImageSourcePixelType>::RealType RealType;

  /**  Type of the measure. */
  typedef typename Superclass::MeasureType         MeasureType;

  /**  Type of the derivative. */
  typedef typename Superclass::DerivativeType      DerivativeType;

  /**  Type of the parameters. */
  typedef typename Superclass::ParametersValueType ParametersValueType;
  typedef typename Superclass::ParametersType      ParametersType;
  typedef Array<unsigned int>                      ParametersMaskType;

  /** Transform and interpolator to pass to the delegate image to image metric. */
  typedef Transform<double, Self::FixedImageDimension, Self::FixedImageDimension>
    TransformType;
  typedef typename TransformType::Pointer 
    TransformTypePointer;

  typedef InterpolateImageFunction<FixedImageType, double>  InterpolatorType;
  typedef typename InterpolatorType::Pointer                InterpolatorTypePointer;

  /** Connect the Fixed Image.  */
  itkSetConstObjectMacro(FixedImage, FixedImageType);

  /** Get the Fixed Image. */
  itkGetConstObjectMacro(FixedImage, FixedImageType);

  /** Connect the Moving Image Source.  */
  virtual void SetMovingImageSource(MovingImageSourceType* source);

  /** Get the Moving Image Source. */
  itkGetObjectMacro(MovingImageSource, MovingImageSourceType);

  /** Overrides SetTransform method of super class. We assume all
      necessary transformations can be handled by the moving image source,
      so this method reports that fact and does not overwrite the transform. */
  void SetTransform(TransformType* transform);

  /** Sets/gets the interpolator to use when accessing values from the
      fixed image. */
  itkSetObjectMacro(Interpolator, InterpolatorType);
  itkGetConstObjectMacro(Interpolator, InterpolatorType);

  /** Get the number of pixels considered in the computation. */
  const unsigned long & GetNumberOfPixelsCounted() const;

  /** Set the region over which the metric will be computed. Forward to
   * the ImageToImageMetric. */
  virtual void SetFixedImageRegion(FixedImageRegionType region);

  /** Get the region over which the metric will be computed */
  virtual const FixedImageRegionType & GetFixedImageRegion();

  /** Set the delegate ImageToImageMetric. */
  virtual void SetDelegateMetric(DelegateMetricType* source);

  /** Get the delegate ImageToImageMetric. */
  itkGetConstObjectMacro( DelegateMetric, DelegateMetricType );

  /** Get the derivative of the cost function. The derivative is estimated
      by forward differences. */
  virtual void GetDerivative(const ParametersType& parameters, DerivativeType& derivative) const;

  /** Get the value of the cost function. The parameters argument should
      contain the values of the active parameters only (in order), not the
      full set of parameters. */
  virtual MeasureType GetValue(const ParametersType& parameters) const;
 
  /** Set active parameters for the moving image Source. The parameters
      argument should contain the values of the active parameters only
      (in order), not the full set of parameters. */
  void SetParameters( const ParametersType & parameters ) const;

  /** Return the number of active parameters required by the 
      ParametricImageSource. */
  virtual unsigned int GetNumberOfParameters(void) const;

  /** Enable/Disable parameters for optimization. */
  void EnableParameter(unsigned int parameterIndex);
  void DisableParameter(unsigned int parameterIndex);
  void         SetParameterEnabled(unsigned int parameterIndex, unsigned int enabled);
  unsigned int GetParameterEnabled(unsigned int parameterIndex);

  /** Initialize the Metric by making sure that all the components
   *  are present and plugged together correctly. */
  virtual void Initialize(void) throw ( ExceptionObject );


protected:
  ImageToParametricImageSourceMetric();
  virtual ~ImageToParametricImageSourceMetric();
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** The fixed image and moving image source. */
  FixedImageConstPointer    m_FixedImage;
  MovingImageSourcePointer  m_MovingImageSource;
  
  /** The delegate image-to-image metric. */
  DelegateMetricTypePointer m_DelegateMetric;

  /** Disable spatial registration by fixing an identity transform. */
  TransformTypePointer      m_Transform;

  /** The interpolator. */
  InterpolatorTypePointer   m_Interpolator;

  /** Mask for parameters. */
  ParametersMaskType        m_ParametersMask;

private:
  ImageToParametricImageSourceMetric(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkImageToParametricImageSourceMetric.txx"
#endif

#endif
