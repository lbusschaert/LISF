module ac_utils
  !! A place for various small utilities used in the
  !! AquaCrop source code

  !! A place for various small utilities.
  use iso_fortran_env, only: int8, &
       int32, &
       real32

  implicit none

  integer, parameter :: sp = real32
  !! double precision real kind

  real(sp), parameter :: ac_zero_threshold = 0.000001_sp

  interface roundc
     module procedure roundc_int8
     module procedure roundc_int32
  end interface roundc


contains


  function roundc_int32(x, mold) result(y)
      !! Returns commercial rounds, following Pascal's banker's rules for rounding
      real(sp), intent(in) :: x
          !! Value to be rounded to an integer
      integer(int32), intent(in) :: mold
          !! Integer determining the kind of the integer result
      integer(int32) :: y
      real(sp) :: x_clipped

      ! Check if x is within the tolerated range for int32 before rounding
      if (x > 2147483647._sp) then
          x_clipped = 2147483647._sp
      elseif (x < -2147483648._sp) then
          x_clipped = -2147483648._sp
      else
          x_clipped = x
      end if

      ! Rounding logic based on x_clipped
      if (abs(x_clipped - floor(x_clipped, kind=int32) - 0.5_sp) < epsilon(0._sp)) then
         if (x_clipped > 0) then
            if (mod(abs(trunc(x_clipped)),2) == 0) then
                y = floor(x_clipped, kind=int32)
            else
                y = ceiling(x_clipped, kind=int32)
            end if
         else
            if (mod(abs(trunc(x_clipped)),2) == 0) then
                y = ceiling(x_clipped, kind=int32)
            else
                y = floor(x_clipped, kind=int32)
            end if
         end if
      else
         ! Standard round for values not ending on 0.5
         y = nint(x_clipped, kind=int32)
      end if

  end function roundc_int32


  function roundc_int8(x, mold) result(y)
    !! Returns commercial rounds, following Pascal's banker's rules for rounding
    real(sp), intent(in) :: x
    !! Value to be rounded to an integer
    integer(int8), intent(in) :: mold
    !! Integer determining the kind of the integer result
    integer(int8) :: y

    if (abs(x - floor(x, kind=int32) - 0.5) < epsilon(0.)) then
       if (x > 0) then
          if (mod(abs(trunc(x)),2) == 0) then
             y = floor(x, kind=int8)
          else
             y = ceiling(x, kind=int8)
          end if
       else
          if (mod(abs(trunc(x)),2) == 0) then
             y = ceiling(x, kind=int8)
          else
             y = floor(x, kind=int8)
          end if
       end if
    else !standard round for values not ending on 0.5
       y = nint(x, kind=int8)
    end if
  end function roundc_int8


  function trunc(x) result(y)
    !! Returns the integer part of x, which is always smaller than (or equal to) x
    !! in absolute value.
    real(sp), intent(in) :: x
    integer(int32) :: y

    if (x > 0) then
       y = floor(x, kind=int32)
    else
       y = ceiling(x, kind=int32)
    end if
  end function trunc

end module ac_utils
