/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkParametricGaussianImageSource.txx,v $
  Language:  C++
  Date:      $Date: 2009-07-12 10:52:52 $
  Version:   $Revision: 1.20 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkParametricGaussianImageSource_txx
#define __itkParametricGaussianImageSource_txx

#include "itkParametricGaussianImageSource.h"

namespace itk
{

template <class TOutputImage>
ParametricGaussianImageSource<TOutputImage>
::ParametricGaussianImageSource()
{
}

template <class TOutputImage>
ParametricGaussianImageSource<TOutputImage>
::~ParametricGaussianImageSource()
{
}

template <class TOutputImage>
void 
ParametricGaussianImageSource<TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
}


template<typename TOutputImage>
void
ParametricGaussianImageSource<TOutputImage>
::SetParameters(const ParametersType& parameters) 
{
  unsigned int i;
  const unsigned int dimensions = itkGetStaticConstMacro(OutputImageDimension);
  ArrayType sigma, mean;
  for (i=0; i<dimensions; i++)
    {
      sigma[i] = parameters[i];
      mean[i]  = parameters[dimensions + i];
    }
  this->SetSigma(sigma);
  this->SetMean(mean);
  this->SetScale(parameters[2*dimensions]);
}


template<typename TOutputImage>
typename ParametricGaussianImageSource<TOutputImage>::ParametersType
ParametricGaussianImageSource<TOutputImage>
::GetParameters() const
{
  ParametersType parameters(GetNumberOfParameters());
  ArrayType sigma = this->GetSigma();
  ArrayType mean  = this->GetMean();
  
  unsigned int i;
  const unsigned int dimensions = itkGetStaticConstMacro(OutputImageDimension);
  for (i=0; i<dimensions; i++)
    {
      parameters[i] = sigma[i];
      parameters[dimensions + i] = mean[i];
    }
  parameters[2*dimensions] = this->GetScale();
  
  return parameters;
}


template<typename TOutputImage>
unsigned int
ParametricGaussianImageSource<TOutputImage>
::GetNumberOfParameters() const
{
  return 2*itkGetStaticConstMacro(OutputImageDimension) + 1;
}

} // end namespace itk

#endif
