
 program getminmax

 double precision, dimension(623091) :: a,b,c

 do i = 1,623091
  read(*,*) a(i),b(i),c(i)
 enddo

 print *,'min max a = ',minval(a),maxval(a)
 print *,'min max b = ',minval(b),maxval(b)
 print *,'min max c = ',minval(c),maxval(c)

 end

