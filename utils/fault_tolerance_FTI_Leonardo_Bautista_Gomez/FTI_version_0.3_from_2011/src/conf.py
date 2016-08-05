#! /usr/bin/python
## * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
##* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
##=======================================================================
##Copyright (C) 2009-2010 Leonardo BAUTISTA GOMEZ
##This program is free software; you can redistribute it and/or modify
##it under the terms of the GNU General Public License (GPL) as published
##of the License, or (at your option) any later version.
##
##This program is distributed in the hope that it will be useful,
##but WITHOUT ANY WARRANTY; without even the implied warranty of
##MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##GNU General Public License for more details.
##
##To read the license please visit http://www.gnu.org/copyleft/gpl.html
##=======================================================================
##
## File          : conf.py
## Created on    : 07 Dec 2010
## Author        : Leonardo BAUTISTA G <leobago@matsulab.is.titech.ac.jp>
##
## Last modified : 14 Jan 2011 (09:22:13 PM)
## Author        : Leonardo BAUTISTA G <leobago@matsulab.is.titech.ac.jp>
## Description   : Python functions for FTI library.
##

import ConfigParser
import os

def getConf(filename):
    config = ConfigParser.RawConfigParser()
    config.read(filename)
    grsz = config.getint("Basic", "Group_size")
    blsz = config.getint("Basic", "Block_size")
    wdsz = config.getint("Basic", "Word_size")
    ndsz = config.getint("Basic", "Node_size")
    mtag = config.getint("Basic", "Mpi_tag")
    fail = config.getint("Basic", "Failure")
    heads = config.getint("Basic", "Heads")
    cdir = config.get("Basic", "Ckpt_dir")
    mdir = config.get("Basic", "Meta_dir")

    conf = [grsz, blsz, wdsz, ndsz, mtag, fail, heads, cdir, mdir]
    return conf


def rmMeta(mfn):
    os.remove(mfn)


def setMeta(mfn, sid, fn, fs, mfs):
    config = ConfigParser.RawConfigParser()
    config.read(mfn)
    config.add_section(str(sid))
    config.set(str(sid),"Ckpt_file_maxs",str(mfs))
    config.set(str(sid),"Ckpt_file_size",str(fs))
    config.set(str(sid),"Ckpt_file_name",fn)
    metadataFile = open(mfn, "wa")
    config.write(metadataFile)
    metadataFile.close()
    return 0


def getTopo(mfn, nid):
    config = ConfigParser.RawConfigParser()
    config.read(mfn)
    hname = config.get("Topology", str(nid))

    return [hname]


def setTopo(mfn, nid, hname):
    config = ConfigParser.RawConfigParser()
    config.read(mfn)
    if not config.has_section("Topology") :
        config.add_section("Topology")
    config.set("Topology",str(nid),str(hname))
    metadataFile = open(mfn, "wa")
    config.write(metadataFile)
    metadataFile.close()
    return 0


def getMeta(filename, sid):
    config = ConfigParser.RawConfigParser()
    config.read(filename)
    fn = config.get(str(sid), "Ckpt_file_name")
    fs = config.getint(str(sid), "Ckpt_file_size")
    mfs = config.getint(str(sid), "Ckpt_file_maxs")

    conf = [fn, fs, mfs]
    return conf

