#!/usr/bin/env python
# -*- coding: utf-8 -*-
#-------------------------------------------------------------
#     GMSH mesh convertion for SPECFEM3D
#
#     by Thomas CURTELIN
#       Centrale Marseille, France, July 2012
#
#     Based on Paul Cristini's equivalent script for 2D meshes
#
# For now it can only handle volumes meshed with hexahedra and boundaries meshed with quadrangles
# (boundary 3D elements must have four nodes on the boundary).
# I tested it with the "homogeneous_halfspace" example provided in the distribution.
#
#--------------------------------------------------------------
#
#     /!\ IMPORTANT REMARKS /!\
#     - Only first order hexahedral 3D-elements are handled by this script
#     - Boundary 2D-elements must thus be first order quadrangles
#     - "xmax", "xmin", "ymax", "ymin", "top" and "bottom" boundaries must be defined as physical surfaces in GMSH
#     - Propagation media "M1", "M2", ... must be defined as physical volumes in GMSH
#
#--------------------------------------------------------------
#--------------------------------------------------------------

#     PACKAGES
####################################################
import sys, string, time
from os.path import splitext, isfile
try:
    from numpy import *
except ImportError:
    print "error: package python-numpy is not installed"


###################################################
#     Save file function (ASCII format)
###################################################

def SauvFicSpecfem(Ng, Ct, Var, Fv):
    # Ng is the name of the file to be written
    # Ct is the number of lines to read
    # Var is the name of the variable containing data to be written
    # Fv is data format (%i for indexes and %f for coordinates)
    savetxt(Ng,(Ct,), fmt='%i')
    fd = open(Ng,'a')
    savetxt(fd, Var, fmt=Fv, delimiter=' ')
    fd.close()
    return

###################################################
#     Read mesh function
###################################################

def OuvreGmsh(Dir,Nom):

      #   Get mesh name

      if splitext(Nom)[-1]=='.msh':
            fic=Nom

      elif splitext(Nom)[-1]=='':
            fic=Nom+'.msh'

      else:
            print 'File extension is not correct'
            print 'script aborted'
            sys.exit()

      #
      # Open the file and get the lines
      ####################################################

      f = file(Dir+fic,'r')
      lignes= f.readlines()
      f.close()

      # Locate information (elements, nodes and physical entities)
      ####################################################
      for ii in range(len(lignes)):
            if lignes[ii]=='$Nodes\n': PosNodes=ii
            if lignes[ii]=='$PhysicalNames\n': PosPhys=ii
            if lignes[ii]=='$Elements\n':
                  PosElem=ii
                  break

      # Element type : second order ONLY
      # 2D elt = 4-node quadrangle : GMSH flag = 3
      # 3D elt = 8-node hexahedron : GMSH flag = 5

      # cf. GMSH documentation available online
      ####################################################

      Ngnod, surfElem, volElem = 4, 3, 5
      len2D, len3D = 4, 8                       #     number of nodes per element according to type

      ###################################################
      # PHYSICAL NAMES
      ###################################################

      #   Get physical surfaces names (borders) and physical volumes names (propagation volumes)

      NbPhysNames = int(string.split(lignes[PosPhys+1])[0])

      # Variable type

      dt = dtype([('dimension',int), ('zone', int), ('name', str, 16)])

      PhysCar=zeros((NbPhysNames,), dtype=dt)

      for Ip in range(NbPhysNames):
            Dim = int(string.split(lignes[PosPhys+2+Ip])[0])    #   2D or 3D
            Zon = int(string.split(lignes[PosPhys+2+Ip])[1])    #   Physical number
            Nam = string.split(lignes[PosPhys+2+Ip])[2][1:-1]   #   Name (xmax, xmin ...)
            PhysCar[Ip] = (Dim, Zon, Nam)                       #   Sorting data
            if Nam == 'xmax':
                  Bord_xmax=Zon
            if Nam == 'xmin':
                  Bord_xmin=Zon
            if Nam == 'ymin':
                  Bord_ymax=Zon
            if Nam == 'ymax':
                  Bord_ymin=Zon
            if Nam == 'bottom':
                  Bord_bottom=Zon
            if Nam == 'top' :
                  Bord_top=Zon

      ###################################################
      print 'Physical Names', PhysCar


      ###################################################
      # GMSH file info
      ####################################################

      Ver=float(string.split(lignes[1])[0])

      File_Type=int(string.split(lignes[1])[1])

      Data_Size=int(string.split(lignes[1])[2])

      ####################################################
      # Nodes
      ####################################################

      NbNodes=int(string.split(lignes[PosNodes+1])[0])    #   Total number of nodes

      print 'Number of nodes: ',NbNodes

      Nodes=zeros((NbNodes,4),dtype=float)                #   Array receiving nodes index and coordinates

      for Ninc in range(NbNodes):
            Nodes[Ninc][0] = int(Ninc+1)
            Nodes[Ninc][1:4] = [float(val) for val in (string.split(lignes[PosNodes+2+Ninc])[1:4])]

      # Save in SPECFEM file format
      ####################################################

      SauvFicSpecfem('nodes_coords_file', NbNodes, Nodes, ['%i','%.9f','%.9f','%.9f'])

      ####################################################
      # Elements
      ####################################################

      NbElements=int(string.split(lignes[PosElem+1])[0])      #   Total number of elements

      #   Initializing arrays

      Elements        = empty((NbElements,len3D+1),dtype=int)   #   3D elements
      Milieu          = empty((NbElements,2),dtype=int)       #   Media index
      Elements3DBord  = empty((NbElements),dtype=int)         #   Volume element next to borders
      Elements2D      = empty((NbElements,len2D),dtype=int)   #   Surface elements (borders)
      #---------------------------------------------------------------------------
      Elements2DBordTop       = empty((NbElements,len2D),dtype=int)
      Elements2DBordBottom    = empty((NbElements,len2D),dtype=int)
      Elements3DBordTop       = zeros((NbElements,len2D+1),dtype=int)
      Elements3DBordBottom    = zeros((NbElements,len2D+1),dtype=int)

      Elements2DBordxmin      = empty((NbElements,len2D),dtype=int)
      Elements2DBordxmax      = empty((NbElements,len2D),dtype=int)
      Elements3DBordxmin      = zeros((NbElements,len2D+1),dtype=int)
      Elements3DBordxmax      = zeros((NbElements,len2D+1),dtype=int)

      Elements2DBordymin      = empty((NbElements,len2D),dtype=int)
      Elements2DBordymax      = empty((NbElements,len2D),dtype=int)
      Elements3DBordymin      = zeros((NbElements,len2D+1),dtype=int)
      Elements3DBordymax      = zeros((NbElements,len2D+1),dtype=int)
      #---------------------------------------------------------------------------

      # Initializing run through elements (surfaces and volumes)

      Ninc2D, Ninc3D = 0, 0

      # Initializing run through boundaries

      Ninc2DBordTop, Ninc2DBordBottom, Ninc2DBordxmax, Ninc2DBordxmin, Ninc2DBordymax, Ninc2DBordymin, = 0, 0, 0, 0, 0, 0

      print 'Number of elements: ', NbElements

      for Ninc in range(NbElements):

            #   Line position

            Pos = PosElem+Ninc+2

            #   Element type position on line

            TypElem = int(string.split(lignes[Pos])[1])

            #   Physical entity number position on line

            ZonP    = int(string.split(lignes[Pos])[3])

            # Initializing material index for Materials_file

            Milieu[Ninc3D]= 1

            #       First case : Surface element

            if TypElem==surfElem:
                  #   Get nodes indexes of the surface element
                  Elements2D[Ninc2D] = [int(val) for val in (string.split(lignes[Pos])[6:])]
                  # Choosing boundary
                  if ZonP==Bord_xmax:
                        Elements2DBordxmax[Ninc2DBordxmax] = Elements2D[Ninc2D]
                        Ninc2DBordxmax+=1
                  if ZonP==Bord_xmin:
                        Elements2DBordxmin[Ninc2DBordxmin] = Elements2D[Ninc2D]
                        Ninc2DBordxmin+=1
                  if ZonP==Bord_ymax:
                        Elements2DBordymax[Ninc2DBordymax] = Elements2D[Ninc2D]
                        Ninc2DBordymax+=1
                  if ZonP==Bord_ymin:
                        Elements2DBordymin[Ninc2DBordymin] = Elements2D[Ninc2D]
                        Ninc2DBordymin+=1
                  if ZonP==Bord_top:
                        Elements2DBordTop[Ninc2DBordTop] = Elements2D[Ninc2D]
                        Ninc2DBordTop+=1
                  if ZonP==Bord_bottom:
                        Elements2DBordBottom[Ninc2DBordBottom] = Elements2D[Ninc2D]
                        Ninc2DBordBottom+=1
                  Ninc2D+=1

                #       Second case : Volume element

            elif TypElem==volElem:
                  Elements[Ninc3D,0] = Ninc3D+1
                  Elements[Ninc3D,1:]= [int(val) for val in (string.split(lignes[Pos])[6:])]
                  Milieu[Ninc3D,0] = Ninc3D+1
                  Milieu[Ninc3D,1] = ZonP-6
                  Ninc3D+=1

            else:
                  print "ERROR : wrong element type flag (3 or 5 only)"

                #       Reduce arrays (exclude zeros elements)

      Elements                =   Elements[:Ninc3D,:]
      Milieu                  =   Milieu[:Ninc3D,:]
      Elements2D              =   Elements2D[:Ninc2D,:]
      Elements2DBordxmin      =   Elements2DBordxmin[:Ninc2DBordxmin,:]
      Elements2DBordxmax      =   Elements2DBordxmax[:Ninc2DBordxmax,:]
      Elements2DBordymin      =   Elements2DBordymin[:Ninc2DBordymin,:]
      Elements2DBordymax      =   Elements2DBordymax[:Ninc2DBordymax,:]
      Elements2DBordTop       =   Elements2DBordTop[:Ninc2DBordTop,:]
      Elements2DBordBottom    =   Elements2DBordBottom[:Ninc2DBordBottom,:]

      # Get nodes from 2D boundary elements
      Elements2DBordFlat=ravel(Elements2D)
      NodesBordC=set(Elements2DBordFlat)
      #-------------------------------------------------------
      NodesBordxmax   = set(ravel(Elements2DBordxmax))
      NodesBordxmin   = set(ravel(Elements2DBordxmin))
      NodesBordymax   = set(ravel(Elements2DBordymax))
      NodesBordymin   = set(ravel(Elements2DBordxmin))
      NodesBordTop    = set(ravel(Elements2DBordTop))
      NodesBordBottom = set(ravel(Elements2DBordBottom))
      #-------------------------------------------------------
      ctBord=0
      ctxmax, ctxmin, ctymax, ctymin, ctt, ctb = 0, 0, 0, 0, 0, 0

      for Ct3D in xrange(Ninc3D):
            #   Test if 3D element contains nodes on boundary
            nodes3DcurrentElement = set(Elements[Ct3D,1:])
            if not set.isdisjoint(nodes3DcurrentElement, NodesBordC):   #   True if there is nodes in common
                #   Choose boundary
                  if not set.isdisjoint(nodes3DcurrentElement, NodesBordxmax):
                                #   Nodes in common between 3D current element and boundary
                        rr = set.intersection(nodes3DcurrentElement, NodesBordxmax)
                        if len(rr) != 4:
                              print "WARNING : wrong 2D boundary element type : ONLY QUADRANGLES"
                              print "Size of wrong intersection :"+str(len(rr))
                              print "Nodes :"
                              print rr
                              sys.exit()
                        else:
                              el = concatenate(([Ct3D+1], list(rr)))
                              Elements3DBordxmax[ctxmax,:] = el
                              ctxmax+=1
              if not set.isdisjoint(nodes3DcurrentElement, NodesBordxmin):
                        rr = set.intersection(nodes3DcurrentElement, NodesBordxmin)
                        if len(rr) != 4:
                              print "WARNING : wrong 2D boundary element type : ONLY QUADRANGLES"
                              print "Size of wrong intersection :"+str(len(rr))
                              print "Nodes :"
                              print rr
                              sys.exit()
                        else:
                              el = concatenate(([Ct3D+1], list(rr)))
                              Elements3DBordxmin[ctxmin,:] = el
                              ctxmin+=1
              if not set.isdisjoint(nodes3DcurrentElement, NodesBordymax):
                        rr = set.intersection(nodes3DcurrentElement, NodesBordymax)
                        if len(rr) != 4:
                              print "WARNING : wrong 2D boundary element type : ONLY QUADRANGLES"
                              print "Size of wrong intersection :"+str(len(rr))
                              print "Nodes :"
                              print rr
                              sys.exit()
                        else:
                              el = concatenate(([Ct3D+1], list(rr)))
                              Elements3DBordymax[ctymax,:] = el
                              ctymax+=1
              if not set.isdisjoint(nodes3DcurrentElement, NodesBordymin):
                        rr = set.intersection(nodes3DcurrentElement, NodesBordymin)
                        if len(rr) != 4:
                              print "WARNING : wrong 2D boundary element type : ONLY QUADRANGLES"
                              print "Size of wrong intersection :"+str(len(rr))
                              print "Nodes :"
                              print rr
                              sys.exit()
                        else:
                              el = concatenate(([Ct3D+1], list(rr)))
                              Elements3DBordymin[ctymin,:] = el
                              ctymin+=1
              if not set.isdisjoint(nodes3DcurrentElement, NodesBordTop):
                        rr = set.intersection(nodes3DcurrentElement, NodesBordTop)
                        if len(rr) != 4:
                              print "WARNING : wrong 2D boundary element type : ONLY QUADRANGLES"
                              print "Size of wrong intersection :"+str(len(rr))
                              print "Nodes :"
                              print rr
                              sys.exit()
                        else:
                              el = concatenate(([Ct3D+1], list(rr)))
                              Elements3DBordTop[ctt,:] = el
                              ctt+=1
              if not set.isdisjoint(nodes3DcurrentElement, NodesBordBottom):
                        rr = set.intersection(nodes3DcurrentElement, NodesBordBottom)
                        if len(rr) != 4:
                              print "WARNING : wrong 2D boundary element type : ONLY QUADRANGLES"
                              print "Size of wrong intersection :"+str(len(rr))
                              print "Nodes :"
                              print rr
                              sys.exit()
                        else:
                              el = concatenate(([Ct3D+1], list(rr)))
                              Elements3DBordBottom[ctb,:] = el
                              ctb+=1
      #       Reducing arrays (exclude zeros elements)

      Elements3DBord=Elements3DBord[:ctBord]
      #----------------------------------------------------------------------
      Elements3DBordTop   = Elements3DBordTop[:ctt,:]
      Elements3DBordxmax    = Elements3DBordxmax[:ctxmax,:]
      Elements3DBordxmin    = Elements3DBordxmin[:ctxmin,:]
      Elements3DBordymax = Elements3DBordymax[:ctymax,:]
      Elements3DBordymin  = Elements3DBordymin[:ctymin,:]
      Elements3DBordBottom   = Elements3DBordBottom[:ctb,:]
      #-----------------------------------------------------------------------
      # Save in SPECFEM file format
      SauvFicSpecfem('mesh_file', Ninc3D, Elements, '%i')
      #
      savetxt('materials_file',Milieu, fmt='%i')
      #
      SauvFicSpecfem('free_surface_file', ctt, Elements3DBordTop, '%i')
      #
      SauvFicSpecfem('absorbing_surface_file_xmax', ctxmax, Elements3DBordxmax, '%i')
      #
      SauvFicSpecfem('absorbing_surface_file_xmin', ctxmin, Elements3DBordxmin, '%i')
      #
      SauvFicSpecfem('absorbing_surface_file_ymax', ctymax, Elements3DBordymax, '%i')
      #
      SauvFicSpecfem('absorbing_surface_file_ymin', ctymin, Elements3DBordymin, '%i')
      #
      SauvFicSpecfem('absorbing_surface_file_bottom', ctb, Elements3DBordBottom, '%i')
      return

if __name__=='__main__':
    set_printoptions(precision=6, threshold=None, edgeitems=None, linewidth=200, suppress=None, nanstr=None, infstr=None)
    #
    Fic = sys.argv[1];                          del sys.argv[1]
    #
    OuvreGmsh('',Fic)

