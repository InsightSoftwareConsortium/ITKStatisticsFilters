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
#ifndef itkKurtosisImageFilter_hxx
#define itkKurtosisImageFilter_hxx
#include "itkKurtosisImageFilter.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkOffset.h"
#include "itkProgressReporter.h"

namespace itk
{
template< typename TInputImage>
KurtosisImageFilter< TInputImage>
::KurtosisImageFilter()
{
    m_Excess = false ;
}

template< typename TInputImage>
void
KurtosisImageFilter< TInputImage>
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
  InputRealType moment4 ;
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
          moment4 = this->GetStandardizedMoment(bit,neighborhoodSize,4);
          value = moment4/(moment2*moment2) ;
          if( this->GetBias() )
          {
            it.Set(value);
          }
          else
          {
            it.Set(1.0/(num-2.0)/(num-3.0)*((num*num-1.0)*value-3*(num-1.0)*(num-1.0))+3.0);
          }
          if( m_Excess )
          {
            it.Set(it.Get()-3);
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
