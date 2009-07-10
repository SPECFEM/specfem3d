! ** Copyright 2004,2007 ENSEIRB, INRIA & CNRS
! **
! ** This file is part of the Scotch software package for static mapping,
! ** graph partitioning and sparse matrix ordering.
! **
! ** This software is governed by the CeCILL-C license under French law
! ** and abiding by the rules of distribution of free software. You can
! ** use, modify and/or redistribute the software under the terms of the
! ** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
! ** URL: "http://www.cecill.info".
! ** 
! ** As a counterpart to the access to the source code and rights to copy,
! ** modify and redistribute granted by the license, users are provided
! ** only with a limited warranty and the software's author, the holder of
! ** the economic rights, and the successive licensors have only limited
! ** liability.
! ** 
! ** In this respect, the user's attention is drawn to the risks associated
! ** with loading, using, modifying and/or developing or reproducing the
! ** software by the user in light of its specific status of free software,
! ** that may mean that it is complicated to manipulate, and that also
! ** therefore means that it is reserved for developers and experienced
! ** professionals having in-depth computer knowledge. Users are therefore
! ** encouraged to load and test the software's suitability as regards
! ** their requirements in conditions enabling the security of their
! ** systems and/or data to be ensured and, more generally, to use and
! ** operate it in the same conditions as regards security.
! ** 
! ** The fact that you are presently reading this means that you have had
! ** knowledge of the CeCILL-C license and that you accept its terms.
! **
! ************************************************************
! **                                                        **
! **   NAME       : ptscotchf.h                             **
! **                                                        **
! **   AUTHOR     : Francois PELLEGRINI                     **
! **                                                        **
! **   FUNCTION   : FORTRAN declaration file for the        **
! **                LibScotch static mapping and sparse     **
! **                matrix block ordering sequential        **
! **                library.                                **
! **                                                        **
! **   DATES      : # Version 3.4  : from : 04 feb 2000     **
! **                                 to     22 oct 2001     **
! **                # Version 4.0  : from : 16 jan 2004     **
! **                                 to     16 jan 2004     **
! **                # Version 5.0  : from : 26 apr 2006     **
! **                                 to     26 apr 2006     **
! **                                                        **
! ************************************************************

! ** Size definitions for the SCOTCH opaque
! ** structures. These structures must be
! ** allocated as arrays of DOUBLEPRECISION
! ** values for proper padding. The dummy
! ** sizes are computed at compile-time by
! ** program "dummysizes".

        INTEGER SCOTCH_ARCHDIM
        INTEGER SCOTCH_DGRAPHDIM
        INTEGER SCOTCH_DORDERDIM
        INTEGER SCOTCH_GEOMDIM
        INTEGER SCOTCH_GRAPHDIM
        INTEGER SCOTCH_MAPDIM
        INTEGER SCOTCH_MESHDIM
        INTEGER SCOTCH_ORDERDIM
        INTEGER SCOTCH_STRATDIM
        PARAMETER (SCOTCH_ARCHDIM   = 5)
        PARAMETER (SCOTCH_DGRAPHDIM = 1)
        PARAMETER (SCOTCH_DORDERDIM = 1)
        PARAMETER (SCOTCH_GEOMDIM   = 2)
        PARAMETER (SCOTCH_GRAPHDIM  = 12)
        PARAMETER (SCOTCH_MAPDIM    = 13)
        PARAMETER (SCOTCH_MESHDIM   = 15)
        PARAMETER (SCOTCH_ORDERDIM  = 12)
        PARAMETER (SCOTCH_STRATDIM  = 1)
