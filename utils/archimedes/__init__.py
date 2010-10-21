#!python
"""Bootstrap merlin installation

If you want to use merlin in your package's setup.py, just include
this directory as a subdirectory alongside setup.py (using
svn:externals), and add this to the top of your setup.py::

    from archimedes import use_merlin
    use_merlin()

If you want to require a specific version of merlin, set a download
mirror, or use an alternate download directory, you can do so by supplying
the appropriate options to ``use_merlin()``.

This file can also be run as a script to install or upgrade merlin.
"""

import sys, os

# Create a new branch of the Merlin project when introducing
# backwards-incompatible changes.  Adjust the following numbers when
# creating a branch.  Branch numbers follow the CVS convention (but
# with no even/odd rules).

branch          = "1"
endOfBranch     = "2a"
default_url     = "http://cheeseshop.python.org/packages/any/m/merlin/"
reqSpec         = "merlin >= %s, < %s" % (branch, endOfBranch)

def use_merlin(package_index=default_url, to_dir=os.curdir, download_delay=15):
    """Automatically find/download merlin and make it available on sys.path

    `package_index` is a URL (ending with a '/') which contains merlin
    eggs.  `to_dir` is the directory where merlin will be downloaded,
    if it is not already available.  If `download_delay` is specified,
    it should be the number of seconds that will be paused before
    initiating a download, should one be required.  If an older
    version of merlin is installed, this routine will print a message
    to ``sys.stderr`` and raise SystemExit in an attempt to abort the
    calling script.
    """

    def bootstrap():
        egg = download_merlin(package_index, to_dir, download_delay)
        sys.path.insert(0, egg)
        import merlin; merlin.bootstrap_install_from = egg

    def rebootstrap():
        del sys.modules['merlin']
        bootstrap()

    try:
        import merlin
        merlin.require(reqSpec)
    except ImportError:
        bootstrap()
    except AttributeError:
        # 'merlin' is pythia's merlin
        del merlin
        rebootstrap()
    except merlin.VersionConflict:
        # 'merlin' is a merlin from a different branch
        del merlin
        rebootstrap()

def download_merlin(package_index=default_url, to_dir=os.curdir, delay=15):
    """Download merlin from a specified location and return its filename

    `package_index` is a URL (ending with a '/') which contains merlin
    eggs. `to_dir` is the directory where the egg will be downloaded.
    `delay` is the number of seconds to pause before an actual
    download attempt.
    """
    import urllib2, urlparse, shutil, re
    from glob import glob
    
    # Avoid repeated downloads
    candidates = []
    candidateMap = {}
    for egg in glob(os.path.join(to_dir, "*.egg")):
        basename = os.path.basename(egg)
        t = basename[:-4].split('-')
        if len(t) == 2:
            name, version = t
            if name == "merlin" and version.startswith(branch + "."):
                candidates.append(basename)
                candidateMap[basename] = egg
    
    if candidates:
        candidates.sort()
        winner = candidates[-1]
        saveto = candidateMap[winner]
    else:
        src = dst = None
        try:
            from distutils import log
            if not sys.stdout.isatty():
                delay = 0
            if delay:
                log.warn("""
---------------------------------------------------------------------------
This script requires merlin version %s to run (even to display
help).  I will attempt to download it for you (from
%s), but
you may need to enable firewall access for this script first.
I will start the download in %d seconds.

(Note: if this machine does not have network access, please obtain the
most recent version of merlin v%s and place it in this directory
before rerunning this script.)
---------------------------------------------------------------------------""",
                    branch, package_index, delay, branch
                ); from time import sleep; sleep(delay)

            # scan the package index
            candidates = []
            candidateMap = {}
            log.warn("Scanning %s", package_index)
            src = urllib2.urlopen(package_index)
            base = src.url     # handle redirects
            data = src.read()
            HREF = re.compile("""href\\s*=\\s*['"]?([^'"> ]+)""", re.I)
            # this is here to fix emacs' cruddy broken syntax highlighting
            for match in HREF.finditer(data):
                link = urlparse.urljoin(base, match.group(1))
                scheme, server, path, parameters, query, fragment = urlparse.urlparse(link)
                basename = urllib2.unquote(path.split('/')[-1])
                if '#' in basename: basename, fragment = basename.split('#',1)
                if basename.endswith('.egg') and '-' in basename:
                    t = basename[:-4].split('-')
                    if len(t) == 2:
                        name, version = t
                        if name == "merlin" and version.startswith(branch + "."):
                            candidates.append(basename)
                            candidateMap[basename] = link

            if not candidates:
                print >>sys.stderr, (
                    "The required version of merlin (v%s) cannot be found.\n"
                    ) % branch
                sys.exit(3)
            
            candidates.sort()
            winner = candidates[-1]
            url = candidateMap[winner]
            saveto = os.path.join(to_dir, winner)
            
            log.warn("Downloading %s", url)
            src = urllib2.urlopen(url)
            # Read/write all in one block, so we don't create a corrupt file
            # if the download is interrupted.
            data = src.read()
            dst = open(saveto,"wb"); dst.write(data)
        finally:
            if src: src.close()
            if dst: dst.close()
    return os.path.realpath(saveto)

def main(argv):
    """Install or upgrade merlin"""

    try:
        import merlin
    except ImportError:
        egg = None
        try:
            egg = download_merlin(delay=0)
            sys.path.insert(0,egg)
            from merlin import main
            return main(['install'] + list(argv) + [egg])   # we're done here
        finally:
            if egg and os.path.exists(egg):
                os.unlink(egg)

    import merlin
    try:
        merlin.require(reqSpec)
    except merlin.VersionConflict:
        from merlin import main
        main(list(argv)+[download_merlin(delay=0)])
        sys.exit(0) # try to force an exit
    else:
        if argv:
            from merlin import main
            main(['install'] + argv)
        else:
            print "Merlin version",branch,"has been installed."



if __name__=='__main__':
    main(sys.argv[1:])
