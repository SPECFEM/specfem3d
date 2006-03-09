#!/usr/bin/env python


from ConfigParser import SafeConfigParser
import pyre.parsing.locators as locators


class HackedConfigParser(SafeConfigParser):

    """Python's SafeConfigParser, hacked to provide file & line info,
    and extended to populate a Pyre Registry."""

    # This class goes to extreme lengths to extract file and line info
    # from Python's ConfigParser module.  It injects proxy 'file' and
    # 'dict' objects into ConfigParser's code.  As the file is parsed,
    # the proxy objects gain control, spying upon the parsing process.
    
    # We should probably write our own parser... *sigh*

    class FileProxy(object):

        def __init__(self):
            self.fp = None
            self.name = "unknown"
            self.lineno = 0

        def readline(self, size=-1):
            line = self.fp.readline()
            if line:
                self.lineno = self.lineno + 1
            return line

        def __getattr__(self, name):
            return getattr(self.fp, name)

    class Section(dict):

        def __init__(self, sectname, node, fp):
            dict.__init__(self)
            self.node = node
            self.fp = fp
        
        def __setitem__(self, key, value):
            locator = locators.file(self.fp.name, self.fp.lineno)
            self.node.setProperty(key, value, locator)
            dict.__setitem__(self, key, value)

    class SectionDict(dict):

        def __init__(self, root):
            dict.__init__(self)
            self.root = root
            self.fp = HackedConfigParser.FileProxy()
        
        def __contains__(self, sectname):
            # Prevent 'ConfigParser' from creating section
            # dictionaries; instead, create our own.
            if not dict.__contains__(self, sectname):
                node = self._getNode(self.root, sectname.split('.'))
                cursect = HackedConfigParser.Section(sectname, node, self.fp)
                self[sectname] = cursect
            return True
        
        def __setitem__(self, key, value):
            dict.__setitem__(self, key, value)

        def _getNode(self, node, path):
            if len(path) == 0:
                return node
            key = path[0].strip()
            return self._getNode(node.getNode(key), path[1:])
    
    def __init__(self, root, defaults=None):
        SafeConfigParser.__init__(self, defaults)
        self._sections = HackedConfigParser.SectionDict(root)

    def _read(self, fp, fpname):
        self._sections.fp.fp = fp
        self._sections.fp.name = fpname
        self._sections.fp.lineno = 0
        SafeConfigParser._read(self, self._sections.fp, fpname)
        self._sections.fp.__init__()


# end of file
