/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkStandardizedMomentImageFilter_hxx
#define itkStandardizedMomentImageFilter_hxx
#include "itkStandardizedMomentImageFilter.h"

#include "itkNeighborhoodInnerProduct.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkOffset.h"
#include "itkProgressReporter.h"
#include "math.h"

namespace itk
{
template< typename TInputImage>
StandardizedMomentImageFilter< TInputImage>
::StandardizedMomentImageFilter()
{
  m_Bias = true;
}

template< typename TInputImage>
typename StandardizedMomentImageFilter< TInputImage>::InputRealType
StandardizedMomentImageFilter< TInputImage>
::GetStandardizedMoment(ConstNeighborhoodIterator< InputImageType > bit, unsigned int neighborhoodSize, unsigned int moment )
{
  InputRealType average = 0.0;
  InputRealType sum = 0.0;
  InputRealType value ;
  InputRealType num = static_cast<InputRealType>(neighborhoodSize);
  for (unsigned int i = 0; i < neighborhoodSize; ++i )
    {
      average += static_cast< InputRealType >( bit.GetPixel(i) );
    }
  average /= num;
  for (unsigned int i = 0; i < neighborhoodSize; ++i )
    {
    value = static_cast< InputRealType >( bit.GetPixel(i) );
    sum += std::pow(value-average,moment);
    }
  return sum/num;
}

} // end namespace itk

#endif
