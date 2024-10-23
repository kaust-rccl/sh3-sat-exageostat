
program fortran_example
    implicit none
    include 'morse_fortran.h'
    integer, parameter:: dp=kind(0.d0) ! double precision
    integer, parameter :: NCPU=2, NGPU=0
    integer, parameter :: N=500, NRHS=1
    double precision, dimension(N*N)    :: A, Acpy
    double precision, dimension(N*NRHS) :: B, X
    double precision :: anorm, bnorm, xnorm, res, eps=1.11022d-16
    integer :: info
    integer :: UPLO=MorseUpper
    logical :: hres


    ! Initialize MORSE with main parameters
    call MORSE_Init(NCPU, NGPU, info)

    ! generate A matrix with random values such that it is spd
    call MORSE_dplgsy( dfloat(N), MorseUpperLower, N, A, N, 51, info )
    Acpy = A

    ! generate RHS
    call MORSE_dplrnt( N, NRHS, B, N, 5673, info )
    X = B

    call MORSE_dpotrf( UPLO, N, A, N, INFO )
    call MORSE_dpotrs( UPLO, N, NRHS, A, N, X, N, info)

    ! compute norms to check the result
    call MORSE_dlange( MorseInfNorm, N, N, Acpy, N, anorm)
    call MORSE_dlange( MorseInfNorm, N, NRHS, B, N, bnorm)
    call MORSE_dlange( MorseInfNorm, N, NRHS, X, N, xnorm)

    ! compute A*X-B, store the result in B
    call MORSE_dgemm( MorseNoTrans, MorseNoTrans, N, NRHS, N, 1.d0, Acpy, N, X, N, -1.d0, B, N, info)
    call MORSE_dlange( MorseInfNorm, N, NRHS, B, N, res)

    ! if hres = 0 then the test succeed
    ! else the test failed
    hres = .TRUE.
    hres = ( res / N / eps / (anorm * xnorm + bnorm ) > 100.0 )
    print *, "   ||Ax-b||       ||A||       ||x||       ||b|| ||Ax-b||/N/eps/(||A||||x||+||b||)"
    if (hres) then
        print *, res, anorm, xnorm, bnorm, res / N / eps / (anorm * xnorm + bnorm ), "FAILURE"
    else
        print *, res, anorm, xnorm, bnorm, res / N / eps / (anorm * xnorm + bnorm), "SUCCESS"
    endif

    ! Finalize MORSE
    call MORSE_Finalize(info)

end program fortran_example
