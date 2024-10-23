#
# Check Example out_of_core
#

set(TESTLIST
    out_of_core
    )

# OOC tests required to create a directory where the data will be written on disk
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ooc)

# define tests
foreach(test ${TESTLIST})
    add_test(example_ooc_${test} ./${test})
endforeach()
