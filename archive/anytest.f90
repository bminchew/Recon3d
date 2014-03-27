 program test_any
            logical l
            l = any((/.true., .true., .true./))
            print *, l
            call section
            contains
              subroutine section
                integer a(2,3), b(2,3), c(3),d,mv(1),z
                a = 1
                b = 1
                b(2,2) = 2
                print *, any(a .eq. b, 1)
                print *, any(a .eq. b, 2)
                if (any(a.eq.1)) print *,'cool' 
 
                c(1) = 2; c(2) = 3; c(3) = 4; d = 3
                mv = minloc(abs(c(2:)-d))
               print *,mv(1)+1
               
            end subroutine section
          end program test_any
