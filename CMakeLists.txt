project(itkilters)
cmake_minimum_required(VERSION 3.5)

find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})

find_package(ITK)
include(${ITK_USE_FILE})
set( ITK_LIBRARIES_WITHIO ${ITK_LIBRARIES} )

set(Filters standarddeviation skewness kurtosis gabor)

foreach(Filter ${Filters})
  SEMMacroBuildCLI(
      NAME ${Filter}
      EXECUTABLE_ONLY
      TARGET_LIBRARIES ${ITK_LIBRARIES_WITHIO}
  )
endforeach()

