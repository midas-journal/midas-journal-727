/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkParametricGaussianImageSource.h,v $
  Language:  C++
  Date:      $Date: 2009-07-12 10:52:52 $
  Version:   $Revision: 1.20 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkParametricGaussianImageSource_h
#define __itkParametricGaussianImageSource_h

#include "itkImageSource.h"
#include "itkFixedArray.h"
#include "itkGaussianImageSource.h"
#include "itkParametricImageSource.h"
#include "itkSize.h"

namespace itk
{

/** \class ParametricGaussianImageSource
 * \brief Adapter class that addes the interface for the ParametricImageSource
 * to the GaussianImageSource.
 *
 * The adapted parameters are the Mean, standard deviation (Sigma), and Scale.
 *
 * This class was developed by Cory Quammen at the Center for Computer
 * Integrated Systems for Microscopy and Manipulation (http:://www.cismm.org)
 * at the University of North Carolina at Chapel Hill and was supported by
 * NIH grant P41-EB002025-25A1.
 *
 * \ingroup DataSources
 */
template <typename TOutputImage>
class ITK_EXPORT ParametricGaussianImageSource : public GaussianImageSource<TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef ParametricGaussianImageSource     Self;
  typedef GaussianImageSource<TOutputImage> Superclass;
  typedef SmartPointer<Self>                Pointer;
  typedef SmartPointer<const Self>          ConstPointer;

  /** Typedef for the output image PixelType. */
  typedef typename Superclass::OutputImagePixelType OutputImagePixelType;

  /** Typedef to describe the output image region type. */
  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

  /** Spacing typedef support.  Spacing holds the size of a pixel.  The
   * spacing is the geometric distance between image samples. */
  typedef typename Superclass::SpacingType SpacingType;

  /** Origin typedef support.  The origin is the geometric coordinates
   * of the index (0,0). */
  typedef typename Superclass::PointType PointType;
  
  /** Direction typedef support.  The direction is the direction
   * cosines of the image. */
  typedef typename Superclass::DirectionType DirectionType;
  
  /** Dimensionality of the output image */
  itkStaticConstMacro(NDimensions, unsigned int, TOutputImage::ImageDimension);

  /** Type used to store gaussian parameters. */
  typedef typename Superclass::ArrayType ArrayType;

  /** Size type matches that used for images */
  typedef typename Superclass::SizeType         SizeType;
  typedef typename Superclass::SizeValueType    SizeValueType;

  /** Run-time type information (and related methods). */
  itkTypeMacro(ParametricGaussianImageSource,GaussianImageSource);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  typedef double                       ParametersValueType;
  typedef Array< ParametersValueType > ParametersType;   

  /** Define methods defined by ParametricImageSource. */
  void           SetParameters(const ParametersType& parameters);
  ParametersType GetParameters() const;
  unsigned int   GetNumberOfParameters() const;
  
protected:
  ParametricGaussianImageSource();
  ~ParametricGaussianImageSource();
  void PrintSelf(std::ostream& os, Indent indent) const;


private:
  ParametricGaussianImageSource(const ParametricGaussianImageSource&); //purposely not implemented
  void operator=(const ParametricGaussianImageSource&); //purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkParametricGaussianImageSource.txx"
#endif

#endif
