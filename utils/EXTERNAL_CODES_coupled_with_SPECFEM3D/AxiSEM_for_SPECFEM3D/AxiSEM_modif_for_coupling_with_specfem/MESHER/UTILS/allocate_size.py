#!/usr/bin/env python

fname_in = 'pdb.f90.bak'
fname_out = 'pdb.f90'
verbose = True
mindims = 2

f_in = open(fname_in, 'r')
f_out = open(fname_out, 'w')

count = 0

for line in f_in:
    f_out.write(line)

    # find allocate statements
    if line.find('allocate') > 0:
        
        # make sure its allocate and not allocated or deallocate
        if line.find('allocated') != -1 or line.find('deallocate') != -1:
            continue

        # make sure its not a comment line
        if line.lstrip()[0] == '!':
            continue

        if verbose:
            print '------------------------------------------------------------------------------------------'
            print line

        varname = line.split('allocate(')[-1]
        varname = varname.split('(')[0]
        
        dims = line.split('allocate('+ varname + '(')[-1]
        dims = dims.rstrip()
        dims = dims[:-2]

        ndims = dims.count(',') + 1
        
        if ndims < mindims:
            if verbose:
                print 'ndims < mindims, continuing'
            continue

        dimsl = []
        for n in range(ndims):
            dimn = dims.split(',')[n]
            if dimn.count(':') == 1:
                dimstart = dimn.split(':')[0]
                dimend = dimn.split(':')[1]
                dimsl.append('(' + dimend + ' - ' + dimstart + ' + 1)')
            else:
                dimsl.append(dimn)

        dimsm = ' * '.join(dimsl)

        if verbose:
            print 'varname: ', varname
            print 'dims:    ', dims
            print 'ndims:   ', ndims
            print 'dimsl:   ', dimsl

            print 'dimsm:   ', dimsm

        whitespace = line.split('allocate')[0]

        out_string = whitespace + 'write(735,\'(a30, i10, f12.5, " MB")\'), &'
        f_out.write(out_string + '\n')

        out_string = whitespace + '   "' + varname + '", &'
        f_out.write(out_string + '\n')

        out_string =  whitespace + '   ' + dimsm + ', &'
        f_out.write(out_string + '\n')

        out_string =  whitespace + '   ' + dimsm + ' * 4. / 1024.**2'
        f_out.write(out_string + '\n\n')

        count += 1


print 'Added %d write statements to %s, output in %s' % (count, fname_in, fname_out)

f_in.close()
f_out.close()
