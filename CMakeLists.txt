cmake_minimum_required(VERSION 3.18)
project(ElectroproductionMomentsInverter VERSION 2.0.0 LANGUAGES CXX)

# Load conventional install-directory names before defining the optional
# checkout-associated install target and its public header destination.
include(GNUInstallDirs)

# ROOT supplies the I/O, special functions, and Minuit2 numerical minimizer used
# by the core physics library.  Imported targets carry ROOT's own include paths
# and transitive link dependencies, avoiding machine-specific flags here.
find_package(ROOT 6.26 REQUIRED COMPONENTS Core RIO Tree MathCore MathMore Minuit2)

add_library(emi_core
  src/Bootstrap.C
  src/Context.C
  src/Evaluation.C
  src/FixedMoments.C
  src/FitRunner.C
  src/Input.C
  src/LeptoMoments.C
  src/Minimizer.C
  src/Model.C
  src/MassModels/MassModel.C
  src/MassModels/PhotoTest.C
  src/PhotoMoments.C
  src/PolarizedModel.C
  src/RandomMoments.C
)
target_include_directories(emi_core
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_compile_features(emi_core PUBLIC cxx_std_17)
target_link_libraries(emi_core
  PUBLIC
    ROOT::Core
    ROOT::RIO
    ROOT::Tree
    ROOT::MathCore
    ROOT::MathMore
    ROOT::Minuit2
)
target_compile_options(emi_core PRIVATE -Wall -Wextra -Wpedantic)

add_executable(emi app/main.C app/RuntimeSettings.C)
target_include_directories(emi PRIVATE app)
target_compile_definitions(emi PRIVATE EMI_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(emi PRIVATE emi_core)

install(TARGETS emi emi_core
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
