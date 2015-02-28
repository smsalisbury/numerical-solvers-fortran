program main
use file_io
implicit none

! DATA DICTIONARY
integer,dimension(:),allocatable          ::  values
real,dimension(:,:),allocatable           ::  matrix
integer,dimension(:),allocatable          ::  forcing
integer,dimension(9)                      ::  ref_values = 100
integer                                   ::  namelist_file

integer                                   ::  i,values_count

! READ NAMELIST
namelist /PARAMETERS/ref_values
namelist_file = file_open('inputs')
read(10,PARAMETERS)

values_count = 0
do i=1,size(ref_values)
  if (ref_values(i) == 100) then
    exit
  end if
  values_count = values_count + 1
end do

allocate(values(values_count))
do i=1,values_count
  values(i) = ref_values(i)
end do

! MAKE MATRIX
allocate(matrix(values_count,values_count))
matrix = transpose(reshape((/ 1.0_wp, 1.0_wp, 1.0_wp, -1.0_wp, 0.0_wp, 1.0_wp, 0.5_wp, 0.0_wp, 0.5_wp /), shape(matrix)))
forcing = 0
forcing(2) = 1
do i=1,size(matrix,1)
  write(*,*)matrix(i, :)
end do
write(*,*)matrix\real(forcing,wp)




close(namelist_file)

end program main