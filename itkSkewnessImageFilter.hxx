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
#ifndef itkSkewnessImageFilter_hxx
#define itkSkewnessImageFilter_hxx
#include "itkSkewnessImageFilter.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkOffset.h"
#include "itkProgressReporter.h"

namespace itk
{
template< typename TInputImage>
SkewnessImageFilter< TInputImage>
::SkewnessImageFilter()
{}

template< typename TInputImage>
void
SkewnessImageFilter< TInputImage>
::ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread,
                       ThreadIdType threadId)
{
  unsigned int i;

  ZeroFluxNeumannBoundaryCondition< InputImageType > nbc;

  ConstNeighborhoodIterator< InputImageType > bit;
  ImageRegionIterator< OutputImageType >      it;

  // Allocate output
  typename OutputImageType::Pointer output = this->GetOutput();
  typename  InputImageType::ConstPointer input  = this->GetInput();

  // Find the data-set boundary "faces"
  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< InputImageType >::FaceListType faceList;
  NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< InputImageType > bC;
  faceList = bC( input, outputRegionForThread, this->GetRadius() );

  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< InputImageType >::FaceListType::iterator fit;

  // support progress methods/callbacks
  ProgressReporter progress( this, threadId, outputRegionForThread.GetNumberOfPixels() );

  InputRealType moment2 ;
  InputRealType moment3 ;
  InputRealType value ;

  // Process each of the boundary faces.  These are N-d regions which border
  // the edge of the buffer.
  for ( fit = faceList.begin(); fit != faceList.end(); ++fit )
    {
    bit = ConstNeighborhoodIterator< InputImageType >(this->GetRadius(),
                                                      input, *fit);
    unsigned int neighborhoodSize = bit.Size();
    InputRealType num = static_cast< InputRealType >( bit.Size() );

    it = ImageRegionIterator< OutputImageType >(output, *fit);
    bit.OverrideBoundaryCondition(&nbc);
    bit.GoToBegin();

    while ( !bit.IsAtEnd() )
      {
      moment2 = this->GetStandardizedMoment(bit,neighborhoodSize,2);
      if( abs(moment2) <= NumericTraits< InputRealType >::min() )
      {
        it.Set(NumericTraits< InputRealType >::ZeroValue());
      }
      else
      {
          moment3 = this->GetStandardizedMoment(bit,neighborhoodSize,3);
          value = moment3/std::pow(moment2,1.5);
          if( this->GetBias() )
          {
            it.Set(value);
          }
          else
          {
            it.Set(std::sqrt((num-1.0)*num)/(num-2.0)*value);
          }
      }
      ++bit;
      ++it;
      progress.CompletedPixel();
      }
    }
}
} // end namespace itk

#endif
