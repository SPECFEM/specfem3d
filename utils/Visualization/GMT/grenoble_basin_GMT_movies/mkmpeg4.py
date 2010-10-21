#!/usr/bin/env python

# Creates an AVI animation from a sequence of images
# The input images can be any format other than TIFF with ZLIB compression

import os, sys, re
from optparse import OptionParser

##############################################################
# Globals
#mencoder = "/home/rsz/pub/linux/mplayer/MPlayer-1.0pre5/bin/mencoder"

##############################################################
# create input options
def CreateParser():
    usage = "%prog [options] image1 ... imagen\n\n"
    usage += "Creates an AVI animation from a sequence of images.\n\n"
    usage += "The input images can be any format other than TIFF with\n"
    usage += "ZLIB compression.\n\n"
    usage += "By default, the MS-MPEG4 v.2 codec is used\n"
    usage += "so it will playback in MS PowerPoint 2000, 2003, \n"
    usage += "MS Media Player 9+, MPlayer and Mac OS X QuickTime Player\n"
    usage += "(if the Mac MS-MPEG4 codec is installed).  DivX output can\n"
    usage += "be used instead of MS-MPEG4v2 with the [-d|--divx] switch.\n\n"
    usage += "Version 0.0\n"
    usage += "Remik.Lastname@noaa.gov where Lastname = Ziemlinski"

    parser = OptionParser(usage=usage)
    parser.add_option("-b", "--bitrate", dest="bitrate",
                      help="bitrate (defaults to 2^23 = 8388608)", type="int", \
                      metavar="BITRATE", default=8388608)
    parser.add_option("-d", "--divx", dest="divx", action="store_true", \
                      help="enable DivX codec instead of MS-MPEG4v2", default=False)
    parser.add_option("-f", "--fps", dest="fps",
                      help="frames per second (defaults to 10)", type="int", \
                      metavar="FPS", default=10)
    parser.add_option("-o", "--out", dest="out",
                      help="output filename (REQUIRED)", metavar="FILE", type="string")
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose", default=False,
                      help="shows progress")
    return parser

##############################################################
def ConvertToSGI(files,verbose):
    """
    Mencoder likes sgi rgb images, so convert if necessary.
    """
    pid = os.getpid()

    # check if SGI already, not likely
    p = os.popen( 'identify -verbose %s 1> /dev/null 2>&1'%(files[0]) )
    info = p.readlines()
    p.close()

    reg = re.compile("\s*Format:\s*SGI.*")
    for i in info:
        if reg.match(i):
            return

    if verbose:
        print 'Converting images to temporary RGB format...'

    for f in files:
        os.system('convert %s %s.%d.sgi 1> /dev/null 2>&1' % (f, f, pid) )

##############################################################
def main():
    parser = CreateParser()
    (options, args) = parser.parse_args()

    out     = options.out
    fps     = options.fps
    divx    = options.divx
    verbose = options.verbose
    files   = args
    bitrate = options.bitrate
    pid     = os.getpid()

    if not out:
        print >>sys.stderr,'ERROR: [-o|--out=] output filename required'
        sys.exit()

    if not files:
        print >>sys.stderr,'ERROR: input files required'
        sys.exit()

    if verbose:
        if divx:
            print 'DivX codec selected'
        else:
            print 'MS-MPEG4v2 codec selected'

    ConvertToSGI(files,verbose)

    if verbose:
        print 'Encoding...pass 1...'

    if divx:
        opts = "vbitrate=%d:mbd=2:keyint=100:v4mv:vqmin=3:vlelim=-4:vcelim=7:lumi_mask=0.07:dark_mask=0.10:naq:vqcomp=0.7:vqblur=0.2:mpeg_quant" % (bitrate)
        cmd = "mencoder -ovc lavc -lavcopts vcodec=mpeg4:vpass=1:%s -mf type=sgi:fps=%d -nosound -o /dev/null mf://\*.%d.sgi" % (opts, fps, pid)

        if not verbose:
            cmd += " 1> /dev/null 2>&1"

        os.system(cmd)

        if verbose:
            print 'Encoding...pass 2...'

# (Emmanuel Chaljub: this step was not working so I just copy/pasted the encoding pass 1)
#        cmd = "mencoder -ovc lavc -lavcopts vcodec=mpeg4:vpass=2:%s -mf type=sgi:fps=%d -nosound -o %s mf://\*.%d.sgi" % (opts, fps, out, pid)

        cmd = "mencoder -ovc lavc -lavcopts vcodec=mpeg4:vpass=1:%s -mf type=sgi:fps=%d -nosound -o %s mf://\*.%d.sgi" % (opts, fps, out, pid)
        if not verbose:
            cmd += " 1> /dev/null 2>&1"

        os.system(cmd)
    else:
        # msmepg4v2
        opts = "vbitrate=%d:mbd=2:keyint=100:vqblur=1.0:cmp=2:subcmp=2:dia=2:mv0:last_pred=3" % (bitrate)
        cmd = "mencoder -ovc lavc -lavcopts vcodec=msmpeg4v2:vpass=1:%s -mf type=sgi:fps=%d -nosound -o /dev/null mf://\*.%d.sgi" % (opts, fps, pid)

        if not verbose:
            cmd += " 1> /dev/null 2>&1"

        os.system(cmd)

        if verbose:
            print 'Encoding...pass 2...'

#        cmd = "mencoder -ovc lavc -lavcopts vcodec=msmpeg4v2:vpass=2:%s -mf type=sgi:fps=%d -nosound -o %s mf://\*.%d.sgi" % (opts, fps, out, pid)
        cmd = "mencoder -ovc lavc -lavcopts vcodec=msmpeg4v2:vpass=1:%s -mf type=sgi:fps=%d -nosound -o %s mf://\*.%d.sgi" % (opts, fps, out, pid)

        if not verbose:
            cmd += " 1> /dev/null 2>&1"

        os.system(cmd)


    # remove temporary files
    os.system("rm *.%d.sgi divx2pass.log" % (pid) )

    if verbose:
        print 'Done'

##############################################################

if __name__=='__main__':
    main()

#======================================================================
#  U.S. Department of Commerce (DOC) Software License for script
#  developed at the Geophysical Fluid Dynamics Laboratory/NOAA.
#
#  1. Scope of License
#
#   Subject to all the terms and conditions of this license, DOC grants
#   USER the royalty-free, nonexclusive, nontransferable, and worldwide
#   rights to reproduce, modify, and distribute this script
#   developed at the Geophysical Fluid Dynamics Laboratory/NOAA, herein
#   referred to as the Product.
#
#  2. Conditions and Limitations of Use
#
#   Warranties.  Neither the U.S. Government, nor any agency or
#     employee thereof, makes any warranties, expressed or implied,
#     with respect to the Product provided under this License,
#     including but not limited to the implied warranties or
#     merchantability and fitness for any particular purpose.
#
#   Liability.  In no event shall the U.S. Government, nor any agency
#     or employee thereof, be liable for any  direct, indirect, or
#     consequential damages flowing from the use of the Product
#     provided under this License.
#
#   Non-Assignment.  Neither this License nor any rights granted
#     hereunder are transferable or assignable without the explicit
#     prior written consent of DOC.
#
#   Names and Logos.  USER shall not substitute its name or logo for the
#     name or logo of DOC, or any of its agencies, in identification of
#     the Product.
#
#   Export of technology.  USER shall comply with all U.S. laws and
#     regulations restricting the export of the Product to other
#     countries.
#
#   Governing Law.  This License shall be governed by the laws of
#     United States as interpreted and applied by the Federal courts
#     in the District of Columbia.
#
#  3. Term of License This License shall remain in effect as long as USER
#     uses the Product in accordance with Paragraphs 1 and 2.
#======================================================================
