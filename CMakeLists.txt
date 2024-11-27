cmake_minimum_required(VERSION 3.0.2)
project(NRS_Gen_Filter)

find_package(catkin REQUIRED COMPONENTS
  roscpp
)

# 다른 package에서 find package하면 내 include directory 역시 
# 자동으로 할당하게 하기 위함
# 만들어 놓은 library 역시 마찬가지
catkin_package(
  INCLUDE_DIRS include
  LIBRARIES NRS_Gen_Filter_Lib
  CATKIN_DEPENDS roscpp
)

# 내 빌드 과정에서 포함해야하는 directory를 의미
include_directories(${catkin_INCLUDE_DIRS} include)

add_library(NRS_Gen_Filter_Lib src/Gen_filter.cpp)
target_link_libraries(NRS_Gen_Filter_Lib ${catkin_LIBRARIES})