MESSAGE(${scsnv_SOURCE_DIR} " " ${scsnv_BINARY_DIR})

add_library(sctxlib
    "pmixture.cpp"
    "pdebug.cpp"
)

link_directories(${scsnv_BINARY_DIR})

target_link_libraries(sctxlib scsnvlib z pthread lzma bz2) #${HTS_LIBRARY} ${BWA_LIBRARY} ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_CXX_HL_LIBRARIES} ${HDF5_C_HL_LIBRARIES} pthread lzma bz2 curl )

add_executable(sctx sctx.cpp)
target_link_libraries(sctx sctxlib)
