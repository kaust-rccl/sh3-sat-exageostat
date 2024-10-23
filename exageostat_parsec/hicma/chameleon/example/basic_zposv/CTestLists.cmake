#
# Check Example basic_zposv
#

set(TESTLIST
    posv_morse_functions
    posv_users_functions
    )

foreach(prec ${RP_CHAMELEON_PRECISIONS})
    string(TOUPPER ${prec} PREC)
    if (CHAMELEON_PREC_${PREC})
        foreach(test ${TESTLIST})
            add_test(example_basic_${prec}${test} ./${prec}${test})
        endforeach()
    endif()
endforeach()