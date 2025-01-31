program mcmc_hastings
    implicit none
    integer, parameter :: n_samples = 10000
    real :: samples(n_samples)
    real :: current, proposal, acceptance_ratio
    integer :: i

    ! Initialize the current sample
    current = 0.0

    ! Run the MCMC algorithm
    do i = 1, n_samples
        ! Generate a proposal from a symmetric proposal distribution (e.g., نرمال)
        proposal = current + random_normal()

        ! Compute the acceptance ratio (پیشنهاد توزیع هدفی)
        acceptance_ratio = target_distribution(proposal) / target_distribution(current)

        ! Accept or reject the proposal
        if (random_number() < acceptance_ratio) then
            current = proposal
        end if

        ! Store the sample
        samples(i) = current
    end do

    ! Output the samples
    open(unit=10, file='samples.txt', status='replace')
    do i = 1, n_samples
        write(10, *) samples(i)
    end do
    close(10)

contains

    ! Define the target distribution (e.g., توزیع نرمال)
    real function target_distribution(x)
        real :: x
        target_distribution = exp(-0.5 * x**2) / sqrt(2.0 * 3.141592653589793)
    end function target_distribution

    ! Function to generate a random number from a normal distribution
    real function random_normal()
        real :: u1, u2, w, mult
        logical :: is_valid

        ! Generate uniformly distributed random numbers
        call random_number(u1)
        call random_number(u2)

        ! Box-Muller transform to generate standard normal random variable
        w = sqrt(-2.0 * log(u1))
        mult = 2.0 * 3.141592653589793 * u2

        random_normal = w * cos(mult)  ! یا sin(mult) برای نمونه دیگر
    end function random_normal

end program mcmc_hastings
