!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program MCMC
    
    implicit none

    ! Number of chains
    integer, parameter :: nchains = 10000
    
    ! Maximum number of iterations
    integer, parameter :: Tmax = 1000

    ! Number of parameters
    integer, parameter :: nparameters = 2
    
    integer :: i, t
    real(8), dimension(nparameters) :: thetamin, thetamax, delta_theta
    real(8), dimension(nparameters, Tmax) :: theta
    real(8), dimension(nparameters) :: thetastar
    real(8) :: fraction
    real(8) :: bivexp
    real(8) :: alpha, u
    real(8) :: chi2total
    real(8) :: RandomReal, NormalRandom
    real(8), dimension(nparameters) :: theta_chain
    real(8) :: chi2total_chain
    real(8) :: random_normal

    open(unit = 11, file = 'chains.txt')

    call random_seed()

    ! BWU

    ! Initiation of blockwise updating

    ! Min of theta1 and theta2
    thetamin = (/ 0.0, 0.0 /)

    ! Max of theta1 and theta2
    thetamax = (/ 10.0, 10.0 /)

    delta_theta = (/ 2.0, 2.0 /)

    do i = 1, nchains

    ! Storage space for samples
    theta = 0

    theta(1, 1) = RandomReal(thetamin(1), thetamax(1))
    theta(2, 1) = RandomReal(thetamin(2), thetamax(2))

    t = 1
    do while (t < Tmax)
        t = t + 1

        ! Generate proposal

        ! thetastar = (/ RandomReal(thetamin(1), thetamax(1)), RandomReal(thetamin(2), thetamax(2)) /)
        
        ! thetastar = (/  NormalRandom(theta(1, t - 1), delta_theta(1)), NormalRandom(theta(2, t - 1), delta_theta(2)) /)
        
        ! thetastar = (/ theta(1, t - 1) + (2.0*RandomReal(0.0, 1.0) - 1.0)*delta_theta(1), theta(2, t - 1) + (2.0*RandomReal(0.0, 1.0) - 1.0)*delta_theta(2) /)

        thetastar = (/ theta(1, t - 1) + random_normal(), theta(2, t - 1) + random_normal() /)

        fraction = bivexp(thetastar(1), thetastar(2))/ &
        bivexp(theta(1, t - 1), theta(2, t - 1))

        alpha = min(1.0, fraction)
        call random_number(u)

        if (u <= alpha) then
            theta(1, t) = thetastar(1)
            theta(2, t) = thetastar(2)

            ! chi2total = -log(bivexp(theta(1, t), theta(2, t)))

            ! write(11,"(4e25.16)") 1.0, chi2total, theta(1, t), theta(2, t)

            theta_chain(1) = theta(1, t)
            theta_chain(2) = theta(2, t)
        else
            theta(1, t) = theta(1, t - 1)
            theta(2, t) = theta(2, t - 1)
        end if

    end do ! while
    
    chi2total_chain = -log(bivexp(theta_chain(1), theta_chain(2)))

    write(11,"(4e25.16)") 1.0, chi2total_chain, theta_chain(1), theta_chain(2)

    end do ! i

    close(11)

end program MCMC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function bivexp(theta1, theta2)

    implicit none

    real(8) :: bivexp
    real(8) :: theta1, theta2

    real(8), parameter :: lambda1 = 0.5d0
    real(8), parameter :: lambda2 = 0.1d0
    real(8), parameter :: lambda = 0.01d0
    real(8), parameter :: maxval = 8.0d0

    ! bivexp = exp(- (lambda1 + lambda)*theta1 - (lambda2 + lambda)*theta2 - lambda*maxval)

    bivexp = exp(- (theta1 - 2.0)**2 - (theta2 - 7.0)**2)

end function bivexp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function RandomReal(a, b)

    implicit none

    real(8) :: RandomReal
    real(8) :: a, b
    real(8) :: r

    call random_seed()

    call random_number(r)

    RandomReal = a + (b - a)*r

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function NormalRandom(mu, sigma)

    implicit none

    real(8) :: NormalRandom
    real(8) :: mu, sigma
    real(8) :: x, y
    real(8) :: r
    real(8) :: PGaussian

    call random_seed()

    do while (.True.)

        call random_number(r)

        x = mu + (2.0*r - 1.0)*sigma
        
        call random_number(y)

        if (y <= PGaussian(x, mu, sigma)) then
            exit
        end if
    
    end do ! while

    NormalRandom = x

end function NormalRandom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
function PGaussian(x, mu, sigma)

    implicit none

    real(8) :: PGaussian
    real(8) :: x, mu, sigma
    real(8), parameter :: pi = 3.1415926535897932384626433832795d0

    PGaussian = (1.0/(sqrt(2.0*pi)*sigma))*exp(-(x - mu)**2/(2.0*sigma**2))

end function PGaussian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Function to generate a random number from a normal distribution
function random_normal()

implicit none

real(8) :: random_normal
real(8) :: u1, u2, w, mult
logical :: is_valid

! Generate uniformly distributed random numbers
call random_number(u1)
call random_number(u2)

! Box-Muller transform to generate standard normal random variable
w = sqrt(-2.0 * log(u1))
mult = 2.0 * 3.141592653589793 * u2

random_normal = w * cos(mult)

end function random_normal