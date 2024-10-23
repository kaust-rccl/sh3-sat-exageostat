cmake_policy(SET CMP0057 NEW)

set( tests_files synthetic_dmle_test.c   synthetic_smle_test.c synthetic_sdmle_test.c
	 synthetic_dmle_locs_rwrapper_test.c synthetic_dmle_rwrapper_test.c real_csv_dmle_test.c
	 real_csv_smle_test.c real_csv_sdmle_test.c)

if(EXAGEOSTAT_USE_NETCDF)
	list(APPEND tests_files real_nc_dmle_1dlocs_test.c real_nc_dmle_2dlocs_test.c real_nc_sdmle_2dlocs_test.c real_nc_sdmle_1dlocs_test.c)
endif()

foreach(test_src ${tests_files})
    get_filename_component(test_exe ${test_src} NAME_WE)


if(EXAGEOSTAT_USE_CUDA)
#message(WARNING "USING CUDA")
	cuda_add_executable(${test_exe} ${test_src} examples.c)
        set_target_properties(${test_exe} PROPERTIES LINKER_LANGUAGE C )
else()
#message(WARNING "NOT USING CUDA")
	add_executable(${test_exe} ${test_src} examples.c)
endif()

    target_link_libraries(${test_exe} exageostat)
    #set_target_properties(test_${test_exe} PROPERTIES OUTPUT_NAME ${test_exe})
    list(APPEND exe_names ${test_exe})
endforeach()


#foreach(test_src ${tests_files})
#    get_filename_component(test_exe ${test_src} NAME_WE)
#    add_executable(${test_exe} ${test_src} examples.c)
#    target_link_libraries(${test_exe} exageostat)
#    #set_target_properties(test_${test_exe} PROPERTIES OUTPUT_NAME ${test_exe})
#    list(APPEND exe_names ${test_exe})
#endforeach()

