
mpif90 -O3 -check nobounds -xHost -ftz -assume buffered_io -assume byterecl -vec-report3 -implicitnone -warn truncated_source -warn argument_checking -warn declarations -warn alignments -warn ignore_loc -warn usage -shared-intel -mcmodel=medium read_absorbing_interfaces.f90 -o ../../bin/xread_absorbing_interfaces


