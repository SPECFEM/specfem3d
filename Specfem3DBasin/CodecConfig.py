#!/usr/bin/env python


from pyre.inventory.odb.Registry import Registry
from pyre.odb.fs.CodecODB import CodecODB
from HackedConfigParser import HackedConfigParser


class CodecConfig(CodecODB):

    def __init__(self):
        CodecODB.__init__(self, encoding='cfg')
        return

    def _decode(self, shelf):
        root = Registry("root")
        parser = HackedConfigParser(root)
        parser.read(shelf.name)
        shelf['inventory'] = root
        shelf._frozen = True
        return


# end of file
