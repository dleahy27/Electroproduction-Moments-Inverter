cmake_minimum_required(VERSION 3.18)
project(ElectroproductionMomentsInverter VERSION 2.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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
  PUBLIC include
  PRIVATE src
)
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

include(GNUInstallDirs)
install(TARGETS emi emi_core
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

enable_testing()
add_test(NAME cli_help COMMAND emi --help)

add_executable(emi_mass_models_test tests/MassModelsTest.C)
target_include_directories(emi_mass_models_test PRIVATE src)
target_link_libraries(emi_mass_models_test PRIVATE emi_core)
target_compile_options(emi_mass_models_test PRIVATE -Wall -Wextra -Wpedantic)
add_test(
  NAME mass_models
  COMMAND emi_mass_models_test
          "${CMAKE_CURRENT_BINARY_DIR}/photo_test_generation.root"
)
