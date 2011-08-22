/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkImageToParametricImageSourceMetric.txx,v $
  Language:  C++
  Date:      $Date: 2009/09/17 20:30:15 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkImageToParametricImageSourceMetric_txx
#define __itkImageToParametricImageSourceMetric_txx

// First, make sure that we include the configuration file.
// This line may be removed once the ThreadSafeTransform gets
// integrated into ITK.
#include "itkConfigure.h"

#include "itkIdentityTransform.h"
#include "itkImageToParametricImageSourceMetric.h"

namespace itk
{

/**
 * Constructor
 */
template <class TFixedImage, class TMovingImageSource> 
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::ImageToParametricImageSourceMetric()
{
  m_FixedImage        = 0; // has to be provided by the user.
  m_MovingImageSource = 0; // has to be provided by the user.
  m_DelegateMetric    = 0; // has to be provided by the user.
  m_Transform         = IdentityTransform<double, Self::FixedImageDimension>::New(); // immutable
  m_Interpolator      = 0; // has to be provided by the user.
  m_ParametersMask    = ParametersMaskType(0);
}


/**
 * Destructor
 */
template <class TFixedImage, class TMovingImageSource> 
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::~ImageToParametricImageSourceMetric()
{

}


template <class TFixedImage, class TMovingImageSource>
void
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::SetMovingImageSource(MovingImageSourceType* source) {
  if (this->m_MovingImageSource != source) {
    this->m_MovingImageSource = source;
    this->Modified();
    
    // Now reinitialize the parameter mask array to match the number
    // of the parameters in the moving image source.
    m_ParametersMask = ParametersMaskType(source->GetNumberOfParameters());
    
    // Initialize to have all enabled parameters.
    for (unsigned int i = 0; i < m_ParametersMask.Size(); i++) {
      m_ParametersMask[i] = 1;
    }
  }
}


template <class TFixedImage, class TMovingImageSource>
void
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::SetTransform(TransformType* transform) {
  itkWarningMacro(<< "Setting the Transform on an "
                  << "ImageToParametrizedImageSourceMetric has no effect. "
                  << "Default IdentityTransform still in use."); 
}


template <class TFixedImage, class TMovingImageSource>
const unsigned long &
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::GetNumberOfPixelsCounted() const {
  if( !m_DelegateMetric )
    {
    itkExceptionMacro(<<"ImageToImageMetric is not present");
    }


  return m_DelegateMetric->GetNumberOfPixelsCounted();
}


template <class TFixedImage, class TMovingImageSource>
void
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::SetFixedImageRegion(FixedImageRegionType region) {
  if( !m_DelegateMetric )
    {
    itkExceptionMacro(<<"ImageToImageMetric is not present");
    }

  m_DelegateMetric->SetFixedImageRegion(region);
}


template <class TFixedImage, class TMovingImageSource>
const typename ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>::FixedImageRegionType & 
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::GetFixedImageRegion() {
  if( !m_DelegateMetric )
    {
    itkExceptionMacro(<<"ImageToImageMetric is not present");
    }  

  return m_DelegateMetric->GetFixedImageRegion();
}


template <class TFixedImage, class TMovingImageSource>
void
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::SetDelegateMetric(DelegateMetricType* source) {
  if (this->m_DelegateMetric != source) {
    this->m_DelegateMetric = source;
    this->Modified();
    
    m_DelegateMetric->SetTransform(m_Transform);
    m_DelegateMetric->SetInterpolator(m_Interpolator);
  }
}


template <class TFixedImage, class TMovingImageSource>
void
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::GetDerivative(const ParametersType& parameters, DerivativeType& derivative) const {
  // Forward differences approximation to the gradient.
  double h = 1e-9;
  double value = GetValue(parameters);
  ParametersType forwardParameters = parameters;

  derivative = DerivativeType(GetNumberOfParameters());

  for (unsigned int i = 0; i < GetNumberOfParameters(); i++) {
    double previousParameterValue = forwardParameters[i];
    forwardParameters[i] += h;
    derivative[i] = (GetValue(forwardParameters) - value) / h;
    forwardParameters[i] = previousParameterValue;
  }

}


template <class TFixedImage, class TMovingImageSource>
typename ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>::MeasureType
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::GetValue(const ParametersType& parameters) const {

  if (!m_MovingImageSource) {
    itkExceptionMacro(<< "Moving image has not been assigned");
  }

  if (!m_DelegateMetric) {
    itkExceptionMacro(<< "Delegate metric has not been assigned");
  }

  if (!m_Interpolator) {
    itkExceptionMacro(<< "Interpolator has not been assigned");
  }

  // Send the parameters to the parametric image source.
  SetParameters(parameters);

  // Now update the parametric image source.
  m_MovingImageSource->GetOutput()->SetRequestedRegionToLargestPossibleRegion();
  m_MovingImageSource->Update();

  m_DelegateMetric->SetFixedImage(m_FixedImage);
  m_DelegateMetric->SetFixedImageRegion(m_FixedImage->
                                        GetLargestPossibleRegion());

  // Have to set the new moving image in the interpolator manually because
  // the delegate image to image metric does this only at initialization.
  MovingImageSourceOutputImagePointerType movingImage = 
    m_MovingImageSource->GetOutput();
  m_Interpolator->SetInputImage(movingImage);

  // Now we can set the moving image in the image to image metric.
  m_DelegateMetric->SetMovingImage(movingImage);

  MeasureType value = m_DelegateMetric->GetValue(parameters);
  return value;
}


template <class TFixedImage, class TMovingImageSource> 
void
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::SetParameters( const ParametersType & parameters ) const
{
  if( !m_MovingImageSource )
    {
    itkExceptionMacro(<<"Moving image source has not been assigned");
    }

  // Iterate through the parameters mask and set only the enabled parameters
  ParametersType allParameters = m_MovingImageSource->GetParameters();
  int enabledIndex = 0;
  for (unsigned int i = 0; i < allParameters.Size(); i++) {
    if (m_ParametersMask[i])
      allParameters[i] = parameters[enabledIndex++];
  }

  m_MovingImageSource->SetParameters( allParameters );
}


template <class TFixedImage, class TMovingImageSource>
unsigned int
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::GetNumberOfParameters(void) const {
  if( !m_MovingImageSource )
    {
    itkExceptionMacro(<<"Moving image source has not been assigned");
    }

  // Count up the enabled parameters
  unsigned int enabledCount = 0;
  for (unsigned int i = 0; i < m_MovingImageSource->GetNumberOfParameters(); i++) {
    if (m_ParametersMask.GetElement(i))
      enabledCount++;
  }
  
  return enabledCount;
}


template <class TFixedImage, class TMovingImageSource> 
void
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::EnableParameter(unsigned int parameterIndex) {
  m_ParametersMask[parameterIndex] = 1;
}


template <class TFixedImage, class TMovingImageSource> 
void
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::DisableParameter(unsigned int parameterIndex) {
  m_ParametersMask[parameterIndex] = 0;
}


template <class TFixedImage, class TMovingImageSource> 
void
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::SetParameterEnabled(unsigned int parameterIndex, unsigned int enabled) {
  m_ParametersMask[parameterIndex] = enabled;
}


template <class TFixedImage, class TMovingImageSource> 
unsigned int
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::GetParameterEnabled(unsigned int parameterIndex) {
  return m_ParametersMask[parameterIndex];
}


template <class TFixedImage, class TMovingImageSource> 
void
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::Initialize(void) throw ( ExceptionObject )
{

  if( !m_MovingImageSource )
    {
    itkExceptionMacro(<<"MovingImageSource is not present");
    }

  if( !m_FixedImage )
    {
    itkExceptionMacro(<<"FixedImage is not present");
    }

  if( !m_DelegateMetric )
    {
    itkExceptionMacro(<<"ImageToImageMetric is not present");
    }

  // If there are any observers on the metric, call them to give the
  // user code a chance to set parameters on the metric
  this->InvokeEvent( InitializeEvent() );
}


template <class TFixedImage, class TMovingImageSource> 
void
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Fixed  Image: " << m_FixedImage.GetPointer() << std::endl;
  os << indent << "Moving Image Source: " << m_MovingImageSource.GetPointer() << std::endl;
  os << indent << "Image to Image Metric: " << m_DelegateMetric.GetPointer() << std::endl;
  os << indent << "Transform: " << m_Transform.GetPointer() << std::endl;
  os << indent << "Interpolator: " << m_Interpolator.GetPointer() << std::endl;
  os << indent << "ParametersMask: " << m_ParametersMask << std::endl;
}


} // end namespace itk


#endif
